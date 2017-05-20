/*
gcc -o brank brank.c -lm
version 2.2
use less memory
pesudocount = min_data/1000
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#define RED "\033[0;32;31m" 
#define NONE "\033[m"

char *USAGE = "Usage: ./brank -i [FILE] -N [RANKNUM]\nOptions:\n-i	Input File\n-N	First N OTU NAME\n-h	Display this information\n";


int GetARGV(int argc, char **argv, char *filename, int *N);
int GetDim(char *filename,int *p,int *n,int N);
void ReadFile(char *filename, double **Mat, char **OTUName, char **GroupInfo,int n);
int RankBetaValue(double **Mat, char **GroupInfo,char **OTUName,int n,int p);
int EnumerateGroupInfo(int *GroupIndex,char **GroupInfo,int n);
double SumAbsBeta(double *array,int num);
double CalculateBetaValue(double **Mat,int i,int j,int *GroupIndex,int GroupNum,int n);
void QuickSort(double *r,int s,int e,char **OTUName);
void NormalizationMatrix(double **Mat,int n,int p);
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
	RankBetaValue(Mat,GroupInfo,OTUName,n,p);

	OutputRankN(OTUName,N);
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

int RankBetaValue(double **Mat, char **GroupInfo,char **OTUName,int n,int p){
	double *beta0;
	beta0=(double*)malloc(p*sizeof(double));
	double *betai;
	betai=(double*)malloc(p*sizeof(double));
	double *ranked_beta;
	ranked_beta=(double*)malloc(p*sizeof(double));
	int *GroupIndex;
	GroupIndex = (int*)malloc(n*sizeof(int));
	int GroupNum = EnumerateGroupInfo(GroupIndex,GroupInfo,n);

	//check GroupNum
	if(GroupNum <= 1) {printf("brank: "RED"fatal error"NONE": the numbers of groups are less than 1\nprogram terminated\n");return -1;}
	NormalizationMatrix(Mat,n,p);

	//Calculate beta of first otu
	beta0[0]=0.0;
	for(int i=1;i<p;i++)
		beta0[i] = CalculateBetaValue(Mat,0,i,GroupIndex,GroupNum,n);
	ranked_beta[0] = SumAbsBeta(beta0,p);
	for(int i=1;i<p;i++){
		for(int j=0;j<p;j++){
			betai[j]=beta0[j]-beta0[i];
		}
		ranked_beta[i] = SumAbsBeta(betai,p);
	}

	QuickSort(ranked_beta,0,p-1,OTUName);
	/*
	for(int i=0;i<p;i++){
		printf("%s\t",OTUName[i]);
		printf("%f\n",ranked_beta[i]);
	}
	*/
	//free mem
	free(beta0);
	free(betai);
	free(ranked_beta);
	free(GroupIndex);
	return 0;
}

double SumAbsBeta(double *array,int num){
	double sum=0;
	for(int i=0;i<num;i++)
		sum+=fabs(array[i]);
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


double CalculateBetaValue(double **Mat,int i,int j,int *GroupIndex,int GroupNum,int n){
	double *SumGroup;
	SumGroup = (double*)malloc(GroupNum*sizeof(double));
	memset(SumGroup,0,GroupNum*sizeof(double));
	int *GroupCount;
	GroupCount = (int*)malloc(GroupNum*sizeof(int));
	memset(GroupCount,0,GroupNum*sizeof(int));
	for(int k=0;k<n;k++){
		SumGroup[GroupIndex[k]] += log(Mat[i][k]/Mat[j][k]);
		GroupCount[GroupIndex[k]] += 1;
	}
	double beta= SumGroup[1]/GroupCount[1] - SumGroup[0]/GroupCount[0];
	//free mem
	free(SumGroup);
	free(GroupCount);
	//TWO groups only
	return beta;
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


void NormalizationMatrix(double **Mat,int n,int p){
	double min_data = 2.0;
	double sum = 0.0;
	for(int j=0;j<n;j++){
		sum = 0.0;
		for(int i=0;i<p;i++)
			sum+=Mat[i][j];
		for(int i=0;i<p;i++){
			Mat[i][j] = Mat[i][j]/sum;
			if(min_data>Mat[i][j]&&Mat[i][j]!=0)
				min_data = Mat[i][j];
		}
	}
	for(int i=0;i<p;i++)
		for(int j=0;j<n;j++)
			if(Mat[i][j] == 0)
				Mat[i][j] = min_data/1000;
}

void OutputRankN(char **OTUName,int N){
	FILE *fp2 = fopen("brank2.txt","w");
	for(int i=0;i<N;i++){
		fwrite(OTUName[i],1,strlen(OTUName[i]),fp2);
		fwrite("\n",1,1,fp2);
	}
	fclose(fp2);

}
