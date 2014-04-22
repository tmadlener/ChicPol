
void random(){
	FILE *SuppFile = fopen("random.txt","w");
	TRandom3 *gRandom = new TRandom3(0);
	for(int i=0; i<68339884; i++){
		fprintf(SuppFile,"%.5f\n", gRandom->Rndm());
	}
	fclose(SuppFile);
}
