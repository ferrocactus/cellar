#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct cmm{
	int cell_index;
	int peak_index;
	float value;
};

int main(int argc, char *argv[]){
	int ncell=atoi(argv[1]);
	float sd=atof(argv[2]);
	FILE *out=fopen(argv[3],"wb");
	int index1;
	int index2;
	int peak;
	int dist;
	int nrow=ncell*(ncell-1)/2;
	float conv[nrow];
	memset(conv,0,sizeof(conv));
	scanf("%d%d%d%d",&index1,&index2,&peak,&dist);
	int previous_peak=peak;
	int index=(index2-1)*ncell+index1-(ncell+ncell-index2+2)*(index2-1)/2;
	conv[index-1]+=pow(0.7788,dist*dist/sd/sd);
	int i;
	struct cmm temp;
	while(scanf("%d%d%d%d",&index1,&index2,&peak,&dist)!=EOF){
		index=(index2-1)*ncell+index1-(ncell+ncell-index2+2)*(index2-1)/2;
		if(previous_peak!=peak){
			for(i=0;i<nrow;i++){
				if(conv[i]>0.01){
					temp.cell_index=i+1;
					temp.peak_index=previous_peak;
					temp.value=conv[i];
					fwrite(&temp,sizeof(cmm),1,out);
					/*printf("%d %d %.2f\n",i+1,previous_peak,conv[i]);*/
				}
			}
			memset(conv,0,sizeof(conv));

		}
		previous_peak=peak;
		conv[index-1]+=pow(0.7788,dist*dist/sd/sd);

	}
	for(i=0;i<nrow;i++){
		if(conv[i]>0.01){
			temp.cell_index=i+1;
			temp.peak_index=peak;
			temp.value=conv[i];
			fwrite(&temp,sizeof(cmm),1,out);
			/*printf("%d %d %.2f\n",i+1,peak,conv[i]);*/
		}
	}
	fclose(out);

}

