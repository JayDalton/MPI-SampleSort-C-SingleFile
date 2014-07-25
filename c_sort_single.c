#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#define LSIZE 24000000  // Lines of File
#define TSIZE 32  // Structur size
#define FOUT "/home/rhelbig/Schreibtisch/PaSys/SampleSortCpp/twitter_files/twitter.out"
#define FLOG "/home/rhelbig/Schreibtisch/PaSys/SampleSortCpp/twitter_files/twitter.log"

FILE* LogFile;
char* MONTHS[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };
char* FILES[] = {
  "/home/rhelbig/Schreibtisch/PaSys/SampleSortCpp/twitter_files/twitter.data.gross.0"
};

int readNumber(char** lptr) {
  char* ptr = *lptr;
  char* line = *lptr;
  while (*ptr != ' ') {
    ptr++;
  }
  *ptr = 0;
  *lptr = ptr + 1;
  return atoi(line);
}

int readMonth(char** lptr) {
  char* ptr = *lptr;
  char* line = *lptr;
  while (*ptr != ' ') {
    ptr++;
  }
  *ptr = 0;
  *lptr = ptr + 1;
  int i, m;
  for (i = 0, m = 1; i < 12; i++, m++) {
    if (strncmp(line, MONTHS[i], 3) == 0) {
      return m;
    }
  }
  fprintf(stderr, "invalid month: %s\n", line);
  exit(3);
}

int countHits(const char* line, const char* key) 
{
  int n = strlen(key);
  int k = strlen(line) - n;
  int i;
  int hits = 0;
  for (i = 0; i < k; i++, line++) {
    if (*line == *key) {
      if (strncmp(line, key, n) == 0) {
        hits++;
      }
    }
  }
  return hits;
}

union {
  unsigned int integer;
  unsigned char byte[4];
} IntConv;

void printTweet(const char* t) 
{
  fprintf(LogFile, "Tweet[");
  int i;
  for (i = 0; i < TSIZE; i++) {
    int k = t[i];
    fprintf(LogFile, "%2x ", k < 0 ? k + 256 : k);
  }
  fprintf(LogFile, "]\n");
}

void printTweetList(char* liste, int anzahl)
{
  int i;
//  for( i = 0; i < anzahl; i++) {
//    printTweet(&liste[i * TSIZE]);
//  }
}

void writeTweet(char* tweet, const int fn, const int ln, const int hits,
    const int month, const int day, char* line) {
  short* ptr1 = (short*) tweet;
  *ptr1 = (short) fn;
  int* ptr2 = (int*) (tweet + 2);
  *ptr2 = ln;
  *(tweet + 6) = (char) hits;
  *(tweet + 7) = (char) month;
  *(tweet + 8) = (char) day;
  int i;
  int n = TSIZE - 9;
  for (i = strlen(line); i < n; i++)
    line[i] = ' '; // padding
  memcpy(tweet + 9, line, n);
       //printTweet(tweet); printf("\n");
}

void saveTweet(char* tweet, const int fn, const int ln, const int hits, const int month, const int day, char* line) 
{
  short* ptr1 = (short*) tweet;
  *ptr1 = (short) fn;
  int* ptr2 = (int*) (tweet + 2);
  *ptr2 = ln;
  *(tweet + 6) = (char) hits;
  *(tweet + 7) = (char) month;
  *(tweet + 8) = (char) day;
  int i;
  int n = TSIZE - 9;
  for (i = strlen(line); i < n; i++) {
    line[i] = ' '; // padding
  }
  memcpy(tweet + 9, line, n);
}

void writeTweets(char* output, const char* fnam, const int lines) 
{
	FILE* file = fopen(fnam, "w");
  if (file == NULL) {
    fprintf(stderr, "open failed: %s\n", fnam);
  }
	int i;
	char* tweet;
	for (i = 0, tweet = output; i < lines; i++, tweet += TSIZE) {
		short* fnp = (short*) tweet;
		int* lnp = (int*) (tweet + 2);
		fprintf(file, "%d %d\n", *fnp, *lnp);
	}
	fclose(file);
}

void loadTweets(char* input, const char* fnam, const int lines, const char* key, int startLine, int endLine, int proc) 
{
  FILE* file = fopen(fnam, "r");
	
  if (file == NULL) {
    fprintf(stderr, "open failed: %s\n", fnam);
  }

  int i;
  int tsize = 32;
  char buffer[1024];
  char* tweet;

  for (i = 0, tweet = input; i < lines; i++) {
    char* line = fgets(buffer, 1024, file);
    if(i >= startLine && i < endLine){
	
    //printf("Reading Line: %i\n", i);
    if (line == NULL) {
      fprintf(stderr, "error reading line %d\n", i);
      exit(2);
    }
    int fn = readNumber(&line);
    int ln = readNumber(&line);
    int month = readMonth(&line);
    int day = readNumber(&line);
    int hits = countHits(line, key);
    saveTweet(tweet , fn, ln, hits, month, day, line);
    tweet += tsize;
    }
  }
  fclose(file);
}

int compTweets(const void* t1, const void* t2)
{
  return memcmp(t2 + 6, t1 + 6, TSIZE - 6);
}

void copyTweet(char *from, char *to) 
{
  memcpy(to, from, TSIZE * sizeof(unsigned char));
}

void copyTweetList(char *from, char *to, int anzahl)
{
  memcpy(to, from, anzahl * TSIZE * sizeof(unsigned char));
}

/**
 * END FORTE
 */

/**** Function Declaration Section ****/

main (int argc, char *argv[])
{
  /* Variable Declarations */

  int             Numprocs, MyRank, Root = 0;
  int             i,j,k, NoElementsToSort;
  int             count, temp;
  char            *Input, *InputData, *Tweets;
  char            FileName[1024];
  unsigned char   *Splitter, *AllSplitter;
  unsigned char   *Buckets, *BucketBuffer, *LocalBucket;
  unsigned char   *OutputBuffer, *Output;
  double          starttime, endtime;
  FILE            *InputFile, *fp;
  MPI_Status      status; 
  
  /**** Initialising ****/
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &Numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

  int FilesSizeAll = sizeof(FILES) / sizeof(FILES[0]);  // Anzahl aller Dateien
  int FilesSizeNode = 1;          // Anzahl Dateien pro Prozess
  int TweetSizeFile = LSIZE / Numprocs;                            // Anzahl Zeilen pro Datei
printf("TSF: %i, FSN: %i\n", TweetSizeFile, FilesSizeNode);
  int TweetSizeNode = TweetSizeFile * FilesSizeNode;    // Anzahl Tweets pro Node
  int TweetSizeAll = TweetSizeNode * Numprocs;          // Anzahl Tweets insgesamt
    
  if (argc != 2) {
    if (MyRank == 0) { 
      printf(" Usage : run size\n");
    }
    MPI_Finalize();
    exit(0);
  }

  if (( TweetSizeAll % Numprocs) != 0){
    if (MyRank == Root) {
      printf("Number of Elements are not divisible by Numprocs \n");
      fprintf(LogFile, "Sending Data: Number of Elements are not divisible by Numprocs \r\n");
    }
    MPI_Finalize();
    exit(0);
  }

  starttime = MPI_Wtime();

  // Init Logging
  sprintf(FileName, "%s.%s.%d", FLOG, argv[1], MyRank);
  LogFile = fopen(FileName, "w");
  
  fprintf(LogFile, "********* Statistics *********\n");
  fprintf(LogFile, "Numprocs: %d\n", Numprocs);
  fprintf(LogFile, "MyRank: %d\n", MyRank);
  fprintf(LogFile, "FilesSizeAll: %i\n", FilesSizeAll);
  fprintf(LogFile, "FilesSizeNode: %i\n", FilesSizeNode);
  fprintf(LogFile, "TweetSizeFile: %i\n", TweetSizeFile);
  fprintf(LogFile, "TweetSizeNode: %i\n", TweetSizeNode);
  fprintf(LogFile, "TweetSizeNode: %i Bytes\n", TweetSizeNode * TSIZE);
  fprintf(LogFile, "TweetSizeAll: %i\n", TweetSizeAll);
  fprintf(LogFile, "TweetSizeAll: %i Bytes\n", TweetSizeAll * TSIZE);
  fprintf(LogFile, "Start Time: %f\n", starttime - MPI_Wtime());
  fprintf(LogFile, "********* Statistics *********\n\n");
    
  /**** Reading Input Parallel ****/
  fprintf(LogFile, "Each Node reading local Input: %i Bytes\n", TweetSizeNode * TSIZE);
  printf("TweetSizeNode: %i\n", TweetSizeNode);
  InputData = (unsigned char*) calloc (TweetSizeNode, TSIZE * sizeof(unsigned char));

  if (InputData == NULL) {
    fprintf(LogFile, "Error : Can not allocate memory \n");
  }

  for (i = 0; i < FilesSizeNode; i++)
  {
    int fileID = MyRank + i * (Numprocs);
    int startLine = MyRank * (LSIZE / Numprocs);
    int endLine = startLine + (LSIZE / Numprocs);
    printf("Rank: %i Startline %i, Endline %i\n", MyRank, startLine, endLine);
   
    fprintf(LogFile, "File: %i ID: %i Name: %s\n", i, fileID, FILES[fileID]);
    loadTweets(InputData + (i * TweetSizeFile * TSIZE), FILES[0], LSIZE, argv[1], startLine, endLine, MyRank);
  }
  
  fprintf(LogFile, "Local Tweet List: %i\n", TweetSizeNode);
  printTweetList(InputData, TweetSizeNode);
//  endtime = MPI_Wtime();
  fprintf(LogFile, "Time Ende File Reading: %f\n", MPI_Wtime() - starttime);
  starttime = MPI_Wtime();
  
  /**** Sorting Locally ****/
  qsort ((unsigned char *)InputData, TweetSizeNode, TSIZE, compTweets);
  fprintf(LogFile, "InputData sorted: %i\n", TweetSizeNode);
  printTweetList(InputData, TweetSizeNode);

  /**** Choosing Local Splitters ****/
  int NumSplitter = Numprocs - 1;
  int TweetBucketSize = TweetSizeNode / Numprocs;

  fprintf(LogFile, "* Choosing Local Splitters: %d, TweetBucketSize: %d  * \n", NumSplitter, TweetBucketSize);
  Splitter = (unsigned char *) calloc (NumSplitter, TSIZE * sizeof(unsigned char));

  for (i = 0; i < NumSplitter; i++) 
  {
    copyTweet(InputData + (i + 1) * (TweetBucketSize * TSIZE), Splitter + (i * TSIZE));
  } 

  fprintf(LogFile, "Splitter local: %d \n",  NumSplitter);
  printTweetList(Splitter, NumSplitter);

  /**** Gathering Local Splitters at Root ****/
  fprintf(LogFile, "* Gathering Local Splitters at Root * \n");
  AllSplitter = (unsigned char*) calloc (Numprocs * NumSplitter, TSIZE * sizeof (unsigned char));

  MPI_Gather (
    Splitter, NumSplitter * TSIZE, MPI_UNSIGNED_CHAR, 
    AllSplitter, NumSplitter * TSIZE, MPI_UNSIGNED_CHAR, 
    Root, MPI_COMM_WORLD
  );

  fprintf(LogFile, "AllSplitter: local %i\n", Numprocs * NumSplitter);
  printTweetList(AllSplitter, Numprocs * NumSplitter);

  /**** Choosing Global Splitters ****/
  if (MyRank == Root) {

    qsort ((unsigned char *) AllSplitter, Numprocs * NumSplitter, TSIZE, compTweets);

    fprintf(LogFile, "All Splitter Root: sorted %i\n", Numprocs * NumSplitter);
    printTweetList(AllSplitter, Numprocs * NumSplitter);

    for (i = 0; i < NumSplitter; i++) 
    {
      copyTweet(AllSplitter + (i + 1) * (NumSplitter * TSIZE), Splitter + (i * TSIZE));
    } 

    fprintf(LogFile, "Splitter Root selected: %i\n", NumSplitter);
    printTweetList(Splitter, NumSplitter);
  }
  
  /**** Broadcasting Global Splitters ****/
  fprintf(LogFile, "Broadcasting Global Splitters\n");
  MPI_Bcast (Splitter, NumSplitter * TSIZE, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  fprintf(LogFile, "Splitter local: %i\n", NumSplitter);
  printTweetList(Splitter, NumSplitter);

  /**** Creating Numprocs Buckets locally ****/
  fprintf(LogFile, "Creating Numprocs Buckets locally: %i + %i\n", TweetSizeNode, Numprocs);
  Buckets = (unsigned char *) calloc (TweetSizeAll + Numprocs, TSIZE * sizeof(unsigned char));

  i = 0;  // Index InputData
  j = 0;  // Index Splitter/Bucket
  k = 1;  // Anzahl Tweets pro Bucket

  for (i = 0; i < TweetSizeNode; i++) {
    
    if (j < NumSplitter) {
      
      if (0 > compTweets(InputData + (i * TSIZE), Splitter + (j * TSIZE))) {

        copyTweet(InputData + (i * TSIZE), Buckets + (((TweetSizeNode + 1) * j) + k) * TSIZE);
        k++;
      
      } else {
        IntConv.integer = k - 1;
        fprintf(LogFile, "Anzahl: %i ", IntConv.integer);
        memcpy(Buckets + ((TweetSizeNode + 1) * j * TSIZE), IntConv.byte, 4);
		    k=1;
			  j++;
		    i--;
       }
    } else {
      copyTweet(InputData + (i * TSIZE), Buckets + (((TweetSizeNode + 1) * j) + k) * TSIZE);
      k++;
	  }
  }
  IntConv.integer = k - 1;
  fprintf(LogFile, "Anzahl: %i\n", IntConv.integer);
  memcpy(Buckets + ((TweetSizeNode + 1) * j * TSIZE), IntConv.byte, 4);

  fprintf(LogFile, "Local Buckets List: %i + %i\n", TweetSizeNode, Numprocs);
  printTweetList(Buckets, TweetSizeAll + Numprocs);

  /**** Sending buckets to respective processors ****/
  fprintf(LogFile, "Sending buckets to respective processors\n");
  BucketBuffer = (unsigned char *) calloc (TweetSizeAll + Numprocs, TSIZE * sizeof (unsigned char));
  
  MPI_Alltoall (  // MPI_Allgather is same?
    Buckets, (TweetSizeNode + 1) * TSIZE, MPI_UNSIGNED_CHAR, 
    BucketBuffer, (TweetSizeNode + 1) * TSIZE, MPI_UNSIGNED_CHAR, 
    MPI_COMM_WORLD
  );
  
  fprintf(LogFile, "BucketBufferAll List: %d\n", TweetSizeAll + Numprocs);
  printTweetList(BucketBuffer, TweetSizeAll + Numprocs);

  /**** Rearranging BucketBuffer ****/
  fprintf(LogFile, "Rearranging BucketBuffer\n");
  LocalBucket = (unsigned char *) calloc (2 * TweetSizeAll,  TSIZE * sizeof(unsigned char));
  fprintf(LogFile, "LocalBucket Alloc: %d\n", 2 * TweetSizeAll);

  count = 1;

  for (j = 0; j < Numprocs; j++) 
  {
    k = 1;
    memcpy(IntConv.byte, BucketBuffer + (TweetSizeNode + 1) * j * TSIZE, 4);
    fprintf(LogFile, "Anzahl: %i ", IntConv.integer);
    for (i = 0; i < IntConv.integer; i++) 
    {
      copyTweet(BucketBuffer + (((TweetSizeNode + 1) * j + k)) * TSIZE, LocalBucket + (count * TSIZE));
      count++;
      k++;
    }
  }
  IntConv.integer = count - 1;
  fprintf(LogFile, "Anzahl: %i\n", IntConv.integer);
  memcpy(LocalBucket, IntConv.byte, 4);

  fprintf(LogFile, "LocalBucket List: %d\n", 2 * TweetSizeAll);
  printTweetList(LocalBucket, 2 * TweetSizeAll);


  /**** Sorting Local Buckets using Bubble Sort ****/
  /*qsort ((char *) InputData, NoofElements_Bloc, sizeof(int), intcompare); */
  fprintf(LogFile, "Sorting Local Buckets using Bubble Sort\n");

  NoElementsToSort = LocalBucket[0];
  memcpy(IntConv.byte, LocalBucket, 4);
  NoElementsToSort = IntConv.integer;
  fprintf(LogFile, "NoElementsToSort: %i\n", NoElementsToSort);

  qsort ((unsigned char *) LocalBucket + TSIZE, NoElementsToSort, TSIZE, compTweets); 
  fprintf(LogFile, "LocalBucket List sorted: %d\n", 2 * TweetBucketSize);

  endtime = MPI_Wtime();
  fprintf(LogFile, "Time Ende Sort Data: %f\n", endtime - starttime);
  starttime = MPI_Wtime();

    /**** Printng the output ****/
  fprintf(LogFile, "**** Printng the output ****\n");
  printTweetList(LocalBucket, 2 * TweetBucketSize);

	sprintf(FileName, "%s.%s.%d", FOUT, argv[1], MyRank);
  writeTweets(LocalBucket + TSIZE, FileName, IntConv.integer);

  endtime = MPI_Wtime();
  fprintf(LogFile, "Time Ende File Writing: %f\n", endtime - starttime);
  starttime = MPI_Wtime();

  free(InputData);
  free(Splitter);
  free(AllSplitter);
  free(Buckets);
  free(BucketBuffer);
  free(LocalBucket);

  fclose(LogFile);

  /**** Finalize ****/
  MPI_Finalize();
}
