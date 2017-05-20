/*
gcc -o frankn frankn.c -lm
version 2.1
use less memory
pesudocount = min_data/1000
if F=inf then F=max_F*2
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#define RED "\033[0;32;31m" 
#define NONE "\033[m"

char *USAGE = "Usage: ./frank -i [FILE] -N [RANKNUM]\nOptions:\n-i	Input File\n-N	First N OTU NAME\n-h	Display this information\n";


int GetARGV(int argc, char **argv, char *filename, int *N);
int GetDim(char *filename,int *p,int *n,int N);
void ReadFile(char *filename, double **Mat, char **OTUName, char **GroupInfo,int n);
int RankFValue(double **Mat, char **GroupInfo,char **OTUName,int n,int p);
double SumFvalue(double *array,int num);
int EnumerateGroupInfo(int *GroupIndex,char **GroupInfo,int n);
double CalculateFValue(double **Mat,int i,int j,int *GroupIndex,int GroupNum,int n);
double average(double *vecs,int num);
void QuickSort(double *r,int s,int e,char **OTUName);
void ReplaceZeroCount(double **Mat,int n,int p);
void OutputRankN(char **OTUName,int N);		

int main(int argc, char **argv){
	char *filename;
	filename = (char *)malloc(200*sizeof(char));
	int N = 0;
	// check whether parameters are too few
	if(GetARGV(argc,argv,filename,&N) == -1) return -1;
	double **Mat;
	char **OTUName;
	char **GroupInfo;
	int p = 0;
	int n = 0;
	//check dimenssions
	if(GetDim(filename,&p,&n,N) == -1) return -1;
	Mat=(double**)malloc(p*sizeof(double*));
	for(int i=0;i<p;i++)
		Mat[i]=(double *)malloc(n*sizeof(double));
	OTUName=(char**)malloc(p*sizeof(char*));
	for(int i=0;i<p;i++)
		OTUName[i]=(char*)malloc(100*sizeof(char));
	GroupInfo=(char**)malloc(n*sizeof(char*));
	for(int i=0;i<n;i++)
		GroupInfo[i]=(char*)malloc(100*sizeof(char));
	ReadFile(filename, Mat, OTUName, GroupInfo, n);
	//check groupnum
	RankFValue(Mat,GroupInfo,OTUName,n,p);

	//OutputRankN(OTUName,N);
	OutputRankN(OTUName,N);

	return 0;
}

int GetARGV(int argc,char **argv,char *filename, int *N){
	int ch;
	int iflag = 0;
	int Nflag =0;
	int hflag = 0;
	while((ch = getopt(argc,argv,"i:N:h"))!= -1 )
	switch(ch){
		case 'i':
			strcpy(filename,optarg);iflag = 1;break;
		case 'N':
			*N = atoi(optarg);Nflag = 1;break;
		case 'h':
			printf("%s\n",USAGE);hflag = 1;break;
		default:
			break;
	}

	//error dealing
	if (hflag == 1){
		return -1;
	}
	else{
		if (iflag == 1 && Nflag == 1) return 0;
		else {printf("brank: "RED"fatal error"NONE": lack of parameters\nprogram terminated.\n");return -1;}
	}

}

int GetDim(char *filename,int *p,int *n,int N){
	FILE *fp = fopen(filename,"r");
	//IO error check
	if(fp == NULL){printf("brank: "RED"fatal error"NONE": IO error\nprogram terminated.\n");return -1;}
	char *line;
	line=(char *)malloc(5000*100*sizeof(char));
	fgets(line,5000*100*sizeof(char),fp);
	for(int i=0;i<strlen(line)-1;i++)
		if(line[i] == '\t' && line[i+1] != '\n')
			*n+=1;
	while(!feof(fp)){
		fgets(line,5000*100*sizeof(char),fp);
		*p+=1;
	}
	*p-=1;
	//free mem
	free(line);
	fclose(fp);

	// error dealing
	if(*p>1 && N<=*p) return 0;
	else{
		if(*p<=1) printf("brank: "RED"fatal error"NONE": the numbers of otus are less than 2\n");
		if(N>*p)  printf("brank: "RED"fatal error"NONE": the parameter N is larger than the number of otus\n"); 
		printf("program terminated.\n");
		return -1;
	}
}

void ReadFile(char *filename, double **Mat, char **OTUName,char **GroupInfo,int n){
	FILE *fp = fopen(filename,"r");
	char *line;
	line=(char*)malloc(n*100*sizeof(char));
	char **p_sample;
	p_sample=(char**)malloc(n*sizeof(char*));
	int col = 0;
	int row = 0;
	char *otuname;
	
	//read first row and get name of group
	fgets(line,n*100*sizeof(char),fp);
	int linelen=strlen(line);
	for(int i=0;i<linelen;i++)
		if(line[i] == '\t' || line[i] == '\n'){
			line[i] = '\0';
			if(i<linelen-1)
				p_sample[col++] = &line[i+1]; 
		}
	for(int i=0;i<n;i++)
		strcpy(GroupInfo[i],p_sample[i]);
	
	//*
	fgets(line,n*100*sizeof(char),fp);
	while(!feof(fp)){
		col = 0;
		linelen=strlen(line);
		otuname = line;
		for(int i=0;i<linelen-1;i++)
			if(line[i] == '\t'){
				line[i] = '\0';
				if(line[i+1]!='\n')
					p_sample[col++] = &line[i+1];
			}
		for(int i=0;i<n;i++)
			Mat[row][i] = atof(p_sample[i]);
		strcpy(OTUName[row],otuname);
		row+=1;
		fgets(line,n*100*sizeof(char),fp);
	}
	//free mem
	free(line);
	free(p_sample);
}

int RankFValue(double **Mat, char **GroupInfo,char **OTUName,int n,int p){
	double *fvalue;
	fvalue=(double*)malloc(p*sizeof(double));
	memset(fvalue,0,p*sizeof(double));
	double *fvalue_sum;
	fvalue_sum = (double*)malloc(p*sizeof(double));
	memset(fvalue_sum,0,p*sizeof(double));

	int GroupIndex[n];
	int GroupNum = EnumerateGroupInfo(GroupIndex,GroupInfo,n);
	ReplaceZeroCount(Mat,n,p);
	//calculate the F value of each pair(p*p)
	double *Max_F_p; //record the max value of F
	Max_F_p = (double*)malloc(p*sizeof(double));
	memset(Max_F_p,0,p*sizeof(double));
	double Max_F = 0.0;
	int *inf_index;
	inf_index=(int*)malloc(p*sizeof(int));
	memset(inf_index,0,p*sizeof(int)); 
	for(int i=0;i<p;i++){
		for(int j=0;j<p;j++){
			if(i!=j){
				fvalue[j] = CalculateFValue(Mat,i,j,GroupIndex,GroupNum,n);
				if(isinf(fvalue[j]) == 0 && Max_F_p[i] <= fvalue[j]) //use max_f*2 instead of inf
					Max_F_p[i] = fvalue[j];
				if(isinf(fvalue[j]) != 0)
					inf_index[i] = 1;
			}
			else
				fvalue[j] = 0.0;
		}
		if(inf_index[i] != 1){
			fvalue_sum[i] = SumFvalue(fvalue,p);
		}
	}
	//find max_f
	for(int i=0;i<p;i++)
		if(Max_F<Max_F_p[i])
			Max_F=Max_F_p[i];
	//calculate otus that have inf F value
	for(int i=0;i<p;i++){
		if(inf_index[i] == 1){
			for(int j=0;j<p;j++){
				if(i!=j){
					fvalue[j]=CalculateFValue(Mat,i,j,GroupIndex,GroupNum,n);
					if(isinf(fvalue[j]) != 0)
						fvalue[j] = Max_F*2;
				}
				else
					fvalue[j]=0.0;

			}
			fvalue_sum[i]=SumFvalue(fvalue,p);
		}
	}
	
	//check GroupNum
	//replace f value (inf & nan)
	QuickSort(fvalue_sum,0,p-1,OTUName);
	free(fvalue);
	free(fvalue_sum);
	free(Max_F_p);
	free(inf_index);
	return 0;
}


double SumFvalue(double *array,int num){
	double sum=0;
	for(int i=0;i<num;i++)
		sum+=array[i];
	return sum;
}

int EnumerateGroupInfo(int *GroupIndex,char **GroupInfo,int n){
	char **p_Group;
	p_Group = (char**)malloc(n*sizeof(char*));
	int GroupNum=1;
	int IndexNum=0;
	GroupIndex[0]=0;
	p_Group[0] = GroupInfo[0];
	for(int i=1;i<n;i++){
		int flag = 0;
		for(int j=0;j<GroupNum;j++){
			if(strcmp(GroupInfo[i],p_Group[j])==0){
				GroupIndex[i] = j;
				flag = 1;
			}
		}
		if(flag == 0){
			GroupIndex[i] = GroupNum;
			p_Group[GroupNum] = GroupInfo[i];
			GroupNum+=1;
		}
	}
	//free mem
	free(p_Group);
	return GroupNum;
}


double CalculateFValue(double **Mat,int i,int j,int *GroupIndex,int GroupNum,int n){
	int GroupCount[GroupNum];
	for(int a=0;a<GroupNum;a++) GroupCount[a] = 0;
	for(int a=0;a<n;a++)
		GroupCount[GroupIndex[a]]+=1;

	double **data;
	data=(double**)malloc(GroupNum*sizeof(double*));
	for(int a=0;a<GroupNum;a++)
		data[a] = (double*)malloc(GroupCount[a]*sizeof(double));
	double **p;
	p=(double**)malloc(GroupNum*sizeof(double*));
	//calculate the log ratio of ith and jth otu
	for(int a=0;a<GroupNum;a++) p[a] = data[a];
	for(int a=0;a<n;a++){
		*p[GroupIndex[a]] = log(Mat[i][a]/Mat[j][a]);
		p[GroupIndex[a]]++;
	}
	//calculate the mean of each group
	double mean[GroupNum];
	for(int a=0;a<GroupNum;a++)
		mean[a] = average(data[a],GroupCount[a]);
	//calculate the mean of all data
	double mean_all = 0;
	int N = 0;
	for(int a=0;a<GroupNum;a++){
		mean_all+=mean[a]*GroupCount[a];
		N+=GroupCount[a];
	}
	mean_all/=N;
	//calculate F value
	double ss_group = 0;
	double ss_error = 0;
	for(int a=0;a<GroupNum;a++)
		ss_group += GroupCount[a]*pow(mean[a]-mean_all,2);
	for(int a=0;a<GroupNum;a++)
		for(int b=0;b<GroupCount[a];b++)
			ss_error+=pow(data[a][b]-mean[a],2);
	double F = (ss_group/(GroupNum-1))/(ss_error/(N-GroupNum));
	//free mem
	for(int a=0;a<GroupNum;a++)
		free(data[a]);
	free(data);
	free(p);
	if(isnan(F)==1)
		return 0.0;
	else
		return F;

}

double average(double *vec,int num){
	double sum=0;
	for(int i=0;i<num;i++)
		sum+=vec[i];
	return sum/num;
}


void QuickSort(double *r,int s,int e,char **OTUName) {
	char *tmpt,*tmpm;
	tmpt = (char*)malloc(100*sizeof(char));
	tmpm = (char*)malloc(100*sizeof(char));
	double t = r[s];
	strcpy(tmpt,OTUName[s]);
	int f = s; 
	int b = e;
	double m = 0; 
	if(s>=e)return;

	while(f<b){ 
		while(r[f]>t) f++;
		while(r[b]<t) b--;
		if(f<=b){
			m = r[f];
			r[f] = r[b];
			r[b] = m;
			strcpy(tmpm,OTUName[f]);
			strcpy(OTUName[f],OTUName[b]);
			strcpy(OTUName[b],tmpm);
			f++;    b--;
		} 
	}
	if(s<b)	QuickSort(r,s,b,OTUName);
	if(f<e) QuickSort(r,f,e,OTUName);
	//free mem
	free(tmpt);
	free(tmpm);
}


void ReplaceZeroCount(double **Mat,int n,int p){
	double min_data = 0.0;
	for(int j=0;j<n;j++){
		min_data = Mat[0][j];
		if(min_data != 0.0)
			break;
	}
	for(int i=0;i<p;i++){
		for(int j=0;j<n;j++){
			if(Mat[i][j] != 0.0 && min_data>Mat[i][j])
				min_data = Mat[i][j];
		}
	}
	for(int i=0;i<p;i++)
		for(int j=0;j<n;j++)
			if(Mat[i][j] == 0)
				Mat[i][j] = min_data/1000;
}

void OutputRankN(char **OTUName,int N){
	FILE *fp2 = fopen("frankn2.txt","w");
	for(int i=0;i<N;i++){
		fwrite(OTUName[i],1,strlen(OTUName[i]),fp2);
		fwrite("\n",1,1,fp2);
	}
	fclose(fp2);

}
