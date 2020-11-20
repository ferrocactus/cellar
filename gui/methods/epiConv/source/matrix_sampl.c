#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct cmm{
	int cell_index;
	int peak_index;
	float value;
};
void retain2uppertri(int ncell, int* retain,int* uppertri){
	int i,j;
	int k=0;int l=1;
	memset(uppertri,0,sizeof(uppertri));
	for(i=2;i<=ncell;i++){
		for(j=1;j<i;j++){
			if(retain[i-1]==1&&retain[j-1]==1){
				uppertri[k]=l;
				l++;
			}
			k++;
		}
	}

}

int main(int argc, char *argv[]){
	FILE *cell_info=fopen(argv[2],"rt");
	char barcode[50];
	char ident[50];
	int ncell=atoi(argv[3]);
	int retain[ncell];
	memset(retain,0,sizeof(retain));
	int uppertri[ncell*(ncell-1)/2];
	memset(uppertri,0,sizeof(uppertri));
	int i=0;
	int j=0;
	while(fscanf(cell_info,"%s%s",barcode,ident)!=EOF){
		if(strcmp(ident,"NA")){
			retain[i]=1;
			j++;
		}
		i++;
	}
	int nvalid=j;
	fclose(cell_info);
	retain2uppertri(ncell,retain,uppertri); /*results stored in static uppertri*/
	/*for(i=0;i<ncell*(ncell-1)/2;i++){
		printf("%d\n",uppertri[i]);
	}*/
	
	FILE *peak_sampl=fopen(argv[4],"rt");
	int npeak,nsampl;
	fscanf(peak_sampl,"%d%d",&npeak,&nsampl);
	int sampl_matrix[npeak][nsampl];
	memset(sampl_matrix,0,sizeof(sampl_matrix));
	int sampl_times[npeak];
	memset(sampl_times,0,sizeof(sampl_times));
	int peak_index,sampl_index;
	while(fscanf(peak_sampl,"%d%d",&peak_index,&sampl_index)!=EOF){
		sampl_matrix[peak_index-1][sampl_times[peak_index-1]]=sampl_index-1;
		sampl_times[peak_index-1]++;
	}
	fclose(peak_sampl);
/*
	for(i=0;i<npeak;i++){
		for(j=0;j<nsampl;j++){
			printf("%d ",sampl_matrix[i][j]);
		}
		printf("\n");
	}
*/

	int nrow_output=nvalid*(nvalid-1)/2;
	float convtmp[nrow_output];
	memset(convtmp,0,sizeof(convtmp));
	int conv_nonzero[nrow_output];
	memset(conv_nonzero,0,sizeof(conv_nonzero));
	int nonzero_count=0;
	float convsum[nrow_output][nsampl];
	memset(convsum,0,sizeof(convsum));
	/*int cell_index;*/
	/*double value;*/
	float peakmean=0;
	float wt;
	int previous_peak;
	FILE *cmat=fopen(argv[1],"rb");
	struct cmm record;
	/*while(fscanf(cmat,"%d%d%lf",&cell_index,&peak_index,&value)!=EOF){*/
	while(fread(&record,sizeof(cmm),1,cmat)!=0){
		if(uppertri[record.cell_index-1]!=0){
			convtmp[uppertri[record.cell_index-1]-1]=record.value;
			conv_nonzero[nonzero_count]=uppertri[record.cell_index-1]-1;
			nonzero_count++;
			peakmean+=record.value;
			previous_peak=record.peak_index;
			break;
		}
	}
	/*while(fscanf(cmat,"%d%d%lf",&cell_index,&peak_index,&value)!=EOF){*/
	while(fread(&record,sizeof(cmm),1,cmat)!=0){
		if(uppertri[record.cell_index-1]==0||sampl_times[record.peak_index-1]==0){
			continue;
		}
		if(previous_peak==record.peak_index){
			convtmp[uppertri[record.cell_index-1]-1]=record.value;
			conv_nonzero[nonzero_count]=uppertri[record.cell_index-1]-1;
			nonzero_count++;
			peakmean+=record.value;
			continue;
		}
		peakmean=peakmean/nrow_output;
		wt=log10(1+1/sqrt(peakmean));
		for(i=0;i<nonzero_count;i++){
			for(j=0;j<sampl_times[previous_peak-1];j++){
				convsum[conv_nonzero[i]][sampl_matrix[previous_peak-1][j]]+=convtmp[conv_nonzero[i]]*wt*wt;
			}
		}
		memset(convtmp,0,sizeof(convtmp));
		memset(conv_nonzero,0,sizeof(conv_nonzero));
		nonzero_count=0;
		peakmean=0;
		previous_peak=record.peak_index;
		convtmp[uppertri[record.cell_index-1]-1]=record.value;
		conv_nonzero[nonzero_count]=uppertri[record.cell_index-1]-1;
		nonzero_count++;
		peakmean+=record.value;

	}
	peakmean=peakmean/nrow_output;
	wt=log10(1+1/sqrt(peakmean));
	for(i=0;i<nonzero_count;i++){
			for(j=0;j<sampl_times[previous_peak-1];j++){
				convsum[conv_nonzero[i]][sampl_matrix[previous_peak-1][j]]+=convtmp[conv_nonzero[i]]*wt*wt;
			}
	}

	for(i=0;i<nrow_output;i++){
		printf("%.4lf",convsum[i][0]);
		for(j=1;j<nsampl;j++){
			printf(" %.4lf",convsum[i][j]);
		}
		printf("\n");
	}


}
