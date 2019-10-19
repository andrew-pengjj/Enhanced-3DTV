/*
 *
 * QccPack: Quantization, compression, and coding utilities
 * Copyright (C) 1997-2009  James E. Fowler
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */
#ifndef MEMC_H
#define MEMC_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>

#include "mex.h"
#define printf mexPrintf

#define TRUE 1
#define FALSE 0

#ifndef MAXINT
#define MAXINT 2147483647
#endif
#ifndef MAXDOUBLE
#define MAXDOUBLE 1.7976931348623157e+308
#endif
#ifndef MAXFLOAT
#define MAXFLOAT 3.40282347e+38F
#endif



#define QCCSTRINGLEN 1200
typedef char QccString[QCCSTRINGLEN + 1];

#define QccFree(p) (p != NULL) ? free(p):(void)0, (p) = NULL

void QccStringCopy(QccString qccstring1, const QccString qccstring2)
{
  if (qccstring1 == NULL)
    return;
  if (qccstring2 == NULL)
    return;

  strncpy((char *)qccstring1, (char *)qccstring2, QCCSTRINGLEN);
  qccstring1[QCCSTRINGLEN] = '\0';
  
}


void QccStringSprintf(QccString qccstring, const char *format, ...)
{
  va_list ap;

  va_start(ap, format);

#ifdef QCC_NO_SNPRINTF
  vsprintf(qccstring, format, ap);
#else
  vsnprintf(qccstring, QCCSTRINGLEN, format, ap);
#endif

  qccstring[QCCSTRINGLEN] = '\0';

  va_end(ap);
}



/*  math.c  */
#define QccMathMax(x, y) (((x) >= (y)) ? (x) : (y))
#define QccMathMin(x, y) (((x) <= (y)) ? (x) : (y))
#define QccMathPercent(x, y) ((y) ? ((double)(x)/(y))*100.0 : (double)0)
#define QccMathModulus(x, y) ((y)?((x)-(y)*((int)floor((double)(x)/(y)))):(x))
#define QccMathLog2(x) ((log(x) / M_LN2))
#define QccMathMedian(a, b, c) (((a) > (b)) ? (((b) > (c)) ? (b) : (((a) > (c)) ? (c) : (a))) : (((a) > (c)) ? (a) : (((b) > (c)) ? (c) : (b))))
#define QCCMATHEPS 3e-16


/*  init.c  */
static QccString QccProgramName;

void QccStringMakeNull(QccString qccstring)
{
  if (qccstring == NULL)
    return;

  qccstring[0] = '\0';
}


void QccConvertToQccString(QccString qccstring, const char *str)
{
  if ((qccstring == NULL) || (str == NULL))
    return;

  strncpy((char *)qccstring, str, QCCSTRINGLEN);
  qccstring[QCCSTRINGLEN] = '\0';
}


void QccExtractProgramName(const char *argv0)
{
  const char *start;

  QccStringMakeNull(QccProgramName);

  if (argv0 == NULL)
    goto Return;

  if (strlen(argv0) < 1)
    goto Return;

  start = strrchr(argv0, '/');
  if (start == NULL)
    QccConvertToQccString((char *)QccProgramName, argv0);
  else
    QccConvertToQccString((char *)QccProgramName, (char *)(start + 1));

 Return:

  return;
}


void QccInit(int argc, char *argv[])
{
  static int first_time = 1;
  QccExtractProgramName((argv != NULL) ? argv[0] : NULL);
  setvbuf(stdout, (char *)NULL, _IOLBF, 0);
}




int QccStringNull(const QccString qccstring)
{
  if (qccstring == NULL)
    return(1);
  return(strlen(qccstring) <= 0);
}


int QccGetProgramName(QccString program_name)
{
  if (program_name == NULL)
    return(0);


  QccStringMakeNull(program_name);
  
  if (QccStringNull(QccProgramName))
    {
      printf("(QccGetProgramName): No valid program name has been set");
      return(1);
    }

  QccStringCopy(program_name, QccProgramName);

  return(0);
  
}



/*  file.c  */
#define QCCFILE_READ 1
#define QCCFILE_WRITE 2
#define QCCFILE_PIPETRUE 1
#define QCCFILE_PIPEFALSE 2


typedef struct
{
  FILE *fileptr;
  int access_mode;
  int pipe_type;
} QccFileEntry;

static int QccFileTableLen = 0;
static QccFileEntry *QccFileTable;





int QccFileExists(const QccString filename)
{
  struct stat statbuf;

  return(!stat(filename, &statbuf));
}



static FILE *QccFileOpenRead(const QccString filename, int *pipe_flag)
{
  int len;
  QccString cmd;
  FILE *infile = NULL;

  *pipe_flag = QCCFILE_PIPEFALSE;
  if (!QccFileExists(filename))
       return(NULL);
 
  if ((infile = fopen(filename, "rb")) == NULL)
      return(NULL);
	len = strlen(filename);
  return(infile);
}
static FILE *QccFileOpenWrite(const QccString filename, int *pipe_flag)
{
  int len;
  QccString cmd;
  FILE *outfile = NULL;

  len = strlen(filename);

  *pipe_flag = QCCFILE_PIPEFALSE;
      if ((outfile = fopen(filename, "wb")) == NULL)
        printf("(QccFileOpenWrite): unable to open %s for writing",
                           filename);
  return(outfile);
}

static int QccFileTableAddEntry(FILE *fileptr,
                                const int access_mode,
                                const int pipe_flag)
{
  int return_value = 0;


  if (QccFileTable == NULL)
    {
      if ((QccFileTable = 
           (QccFileEntry *)malloc(sizeof(QccFileEntry))) == NULL)
        {
          printf("(QccFileTableAddEntry): Error allocating memory");
          goto Error;
        }
      QccFileTableLen = 1;
    }
  else
    {
      QccFileTableLen++;
      if ((QccFileTable = 
           (QccFileEntry *)realloc((void *)QccFileTable,
                                   sizeof(QccFileEntry)*QccFileTableLen)) ==
          NULL)
        {
          printf("(QccFileTableAddEntry): Error reallocating memory");
          goto Error;
        }
    }

  QccFileTable[QccFileTableLen - 1].fileptr = fileptr;
  QccFileTable[QccFileTableLen - 1].access_mode = access_mode;
  QccFileTable[QccFileTableLen - 1].pipe_type = pipe_flag;

  return_value = 0;
  goto Return;

 Error:
  return_value = 1;

 Return:

  return(return_value);
}


FILE *QccFileOpen(const QccString filename, const QccString mode)
{
  FILE *fileptr = NULL;
  int pipe_flag;
  int access_mode;

  if (mode == NULL)
    {
     printf("(QccFileOpen): Error invalid open mode");
      return(NULL);
    }

  if (mode[0] == 'r')
    {
      if ((fileptr = QccFileOpenRead(filename, &pipe_flag)) == NULL)
        {
         printf("(QccFileOpen): Error calling QccFileOpenRead()");
          return(NULL);
        }

      access_mode = QCCFILE_READ;
    }
  else
    if (mode[0] == 'w')
      {
        if ((fileptr = QccFileOpenWrite(filename, &pipe_flag)) == NULL)
          {
            printf("(QccFileOpen): Error calling QccFileOpenWrite()");
            return(NULL);
          }

        access_mode = QCCFILE_WRITE;
      }
    else
      {
        printf("(QccFileOpen): Open mode %s is invalid",
                           mode);
        return(NULL);
      }

  if (QccFileTableAddEntry(fileptr, access_mode, pipe_flag))
    {
     printf("(QccFileOpen): Error calling QccFileTableAddEntry()");
      return(NULL);
    }

  return(fileptr);
}


int QccFileReadString(FILE *infile, QccString s)
{
  QccString cmd;

  if (infile == NULL)
    return(0);

  if (s == NULL)
    return(0);

  QccStringSprintf(cmd, "%%%ds", QCCSTRINGLEN);
  fscanf(infile, cmd, s);

  return(ferror(infile) || feof(infile));
}


static int QccFileSkipComment(FILE *infile)
{
  unsigned char ch;

  if (infile == NULL)
    return(0);

  do
    {
      ch = fgetc(infile);
      if (ferror(infile))
        return(1);
    }
  while (ch != '\n');

  ungetc(ch, infile);

  return(0);
}


int QccFileSkipWhiteSpace(FILE *infile, int skip_comments_flag)
{
  unsigned char ch;

  if (infile == NULL)
    return(0);

  do
    {
      ch = fgetc(infile);
      if (ferror(infile))
        return(1);

      if (((ch == '#') || (ch == '*') || (ch == '/'))
          && skip_comments_flag)
        {
          if (QccFileSkipComment(infile))
            return(1);
          ch = fgetc(infile);
          if (ferror(infile))
            return(1);
        }
    }
  while ((ch == ' ') || (ch == '\n') || (ch == '\t') || (ch == '\r'));

  ungetc(ch, infile);
  return(0);

}

int QccFileReadMagicNumber(FILE *infile, QccString magic_num,
                           int *major_version_number,
                           int *minor_version_number)
{
  QccString header_value;
  int major, minor;

  if (infile == NULL)
    return(0);

  if (magic_num == NULL)
    return(1);

  if (QccFileReadString(infile, header_value))
    return(1);
  if (QccFileSkipWhiteSpace(infile, 1))
    return(1);

  sscanf(header_value, "%[^0-9]%d.%d", magic_num, &major, &minor);

  if (major_version_number != NULL)
    *major_version_number = major;
  if (minor_version_number != NULL)
    *minor_version_number = minor;

  return(0);
}

int QccFileWriteString(FILE *outfile, const QccString s)
{
  if (outfile == NULL)
    return(0);

  if (s == NULL)
    return(0);

  fprintf(outfile, "%s", s);

  return(ferror(outfile));
}


int QccFileWriteMagicNumber(FILE *outfile, const QccString magic_num)
{
  QccString header_value;
  int major;
  int minor;

  if (outfile == NULL)
    return(0);

  QccStringSprintf(header_value, "%s\n",
          magic_num);

  if (QccFileWriteString(outfile, header_value))
    {
      printf("(QccFileWriteMagicNumber): Error writing magic number");
      goto QccErr;
    }

  return(0);

 QccErr:
  return(1);
}

int QccFileGetMagicNumber(const QccString filename, QccString magic_num)
{
  FILE *infile;

  if (filename == NULL)
    return(0);
  if (magic_num == NULL)
    return (0);

  if ((infile = QccFileOpen(filename, "r")) == NULL)
    {
      printf("(QccFileGetMagicNumber): Error opening %s for reading",
                         filename);
      return(1);
    }
  
  if (QccFileReadMagicNumber(infile, magic_num,
                             NULL, NULL))
    {
      printf("(QccFileGetMagicNumber): Error reading magic number of %s",
                         filename);
      return(1);
    }
  
  QccFileClose(infile);
  
  return(0);
}


static int QccFileGetPipeType(FILE *fileptr)
{
  int file_cnt;
  int return_value = 0;

  for (file_cnt = 0; file_cnt < QccFileTableLen; file_cnt++)
    if (QccFileTable[file_cnt].fileptr == fileptr)
      {
        return_value = QccFileTable[file_cnt].pipe_type;
        break;
      }


  return(return_value);
}



static int QccFileTableRemoveEntry(FILE *fileptr)
{
  int file_cnt, file_cnt2, file_cnt3;
  QccFileEntry *tmp_table;
  int return_value = 0;

  for (file_cnt = 0; file_cnt < QccFileTableLen; file_cnt++)
    if (QccFileTable[file_cnt].fileptr == fileptr)
      break;

  if (file_cnt >= QccFileTableLen)
    {
      printf("(QccFileTableRemoveEntry): fileptr is not in table");
      goto Error;
    }

  if (QccFileTableLen == 1)
    {
      QccFree(QccFileTable);
      QccFileTableLen = 0;
      QccFileTable = NULL;
      goto Return;
    }

  if ((tmp_table = 
       (QccFileEntry *)malloc(sizeof(QccFileEntry)*(QccFileTableLen - 1))) == 
      NULL)
    {
      printf("(QccFileTableRemoveEntry): Error allocating memory");
      goto Error;
    }

  for (file_cnt2 = 0, file_cnt3 = 0; file_cnt2 < QccFileTableLen;
       file_cnt2++)
    if (file_cnt2 != file_cnt)
      tmp_table[file_cnt3++] = QccFileTable[file_cnt2];

  QccFree(QccFileTable);
  QccFileTable = tmp_table;
  QccFileTableLen--;

  return_value = 0;
  goto Return;

 Error:
  return_value = 1;

 Return:

  return(return_value);
}


int QccFileGetExtension(const QccString filename, QccString extension)
{
  char *ptr;
  QccString filename2;

  QccStringCopy(filename2, filename);

  if ((ptr = strrchr(filename2, '.')) == NULL)
    {
      QccStringMakeNull(extension);
      return(1);
    }

  if (!strcmp(ptr + 1, "gz"))
    {
      *ptr = '\0';

      if ((ptr = strrchr(filename2, '.')) == NULL)
        {
          QccStringMakeNull(extension);
          return(1);
        }
    }

  QccStringCopy(extension, ptr + 1);

  return(0);
}


int QccFileClose(FILE *fileptr)
{
  int pipe_flag;

  if (fileptr == NULL)
    return(0);

  if (!(pipe_flag = QccFileGetPipeType(fileptr)))
    {
      printf("(QccFileClose): fileptr is not in QccFileTable");
      return(1);
    }

  if (pipe_flag == QCCFILE_PIPETRUE)
    pclose(fileptr);
  else
    fclose(fileptr);

  if (QccFileTableRemoveEntry(fileptr))
    {
      printf("(QccFileClose): Error calling QccFileTableRemoveEntry()");
      return(1);
    }

  return(0);
}

/*  binary_values.c  */
#define QCC_INT_SIZE 4

int QccBinaryCharToInt(const unsigned char *ch, unsigned int *val)
{
  int byte_cnt;

  if ((val == NULL) || (ch == NULL))
    return(0);

  /*  MSB first  */
  *val = 0;
  for (byte_cnt = 0; byte_cnt < QCC_INT_SIZE; byte_cnt++)
    {
      *val <<= 8;
      *val |= ch[byte_cnt];
    }

  return(0);
}


int QccBinaryCharToFloat(const unsigned char *ch, float *val)
{
  /*  Relies on floats being QCC_INT_SIZE bytes long  */
  union
  {
    float d;
    unsigned int i;
  } tmp;

  if ((val == NULL) || (ch == NULL))
    return(0);

  QccBinaryCharToInt(ch, &(tmp.i));
  *val = tmp.d;

  return(0);
}


int QccBinaryIntToChar(unsigned int val, unsigned char *ch)
{
  int byte_cnt;

  if (ch == NULL)
    return(0);

  /*  MSB first  */
  for (byte_cnt = QCC_INT_SIZE - 1; byte_cnt >= 0; byte_cnt--)
    {
      ch[byte_cnt] = (val & 0xff);
      val >>= 8;
    }

  return(0);
}


int QccBinaryFloatToChar(float val, unsigned char *ch)
{
  /*  Relies on floats being QCC_INT_SIZE bytes long  */
  union
  {
    float d;
    unsigned int i;
  } tmp;

  if (ch == NULL)
    return(0);

  tmp.d = val;
  QccBinaryIntToChar(tmp.i, ch);

  return(0);
}


int QccFileReadDouble(FILE *infile, double *val)
{
  unsigned char ch[QCC_INT_SIZE];
  float tmp;

  if (infile == NULL)
    return(0);

  fread(ch, 1, QCC_INT_SIZE, infile);
  if (ferror(infile) || feof(infile))
    return(1);

  if (val != NULL)
    {
      QccBinaryCharToFloat(ch, &tmp);
      
      if (val != NULL)
        *val = tmp;
    }

  return(0);
}

int QccFileWriteDouble(FILE *outfile, double val)
{
  unsigned char ch[QCC_INT_SIZE];
  float tmp;

  if (outfile == NULL)
    return(0);

  tmp = val;
  QccBinaryFloatToChar(tmp, ch);
  fwrite(ch, 1, QCC_INT_SIZE, outfile);
  return(ferror(outfile));

}


/*  vector.c  */
#define QCCVECTOR_SORTASCENDING 1
#define QCCVECTOR_SORTDESCENDING 2
#define QCCVECTOR_EVEN 0
#define QCCVECTOR_ODD 1
typedef double *QccVector;


QccVector QccVectorAlloc(int vector_dimension)
{
  QccVector vector;

  if (vector_dimension <= 0)
    return(NULL);

  if ((vector = (QccVector)calloc(vector_dimension, sizeof(double))) == NULL)
    printf("(QccVectorAlloc): Error allocating memory");

  return(vector);
}

void QccVectorFree(QccVector vector)
{
  if (vector != NULL)
    QccFree(vector);
}


int QccVectorZero(QccVector vector, int vector_dimension)
{
  int component;

  if ((vector == NULL) || (vector_dimension <= 0))
    return(0);
  for (component = 0; component < vector_dimension; component++)
    vector[component] = 0;

  return(0);
}


/*  matrix.c  */
typedef QccVector *QccMatrix;

QccMatrix QccMatrixAlloc(int num_rows, int num_cols)
{
  QccMatrix matrix;
  int row;
  
  if ((num_rows <= 0) || (num_cols <= 0))
    return(NULL);

  if ((matrix = 
       (QccMatrix)malloc(sizeof(QccVector)*num_rows)) == NULL)
    {
      printf("(QccMatrixAlloc): Error allocating memory");
      return(NULL);
    }
  
  for (row = 0; row < num_rows; row++)
    if ((matrix[row] = QccVectorAlloc(num_cols)) == NULL)
      {
        printf("(QccMatrixAlloc): Error allocating memory");
        QccFree(matrix);
        return(NULL);
      }
  
  return(matrix);
}


void QccMatrixFree(QccMatrix matrix, int num_rows)
{
  int row;
  
  if (matrix != NULL)
    {
      for (row = 0; row < num_rows; row++)
        QccVectorFree(matrix[row]);

      QccFree(matrix);
    }
}




/*  filter.c  */
#define QCCFILTER_CAUSAL 0
#define QCCFILTER_ANTICAUSAL 1
#define QCCFILTER_SYMMETRICWHOLE 2
#define QCCFILTER_SYMMETRICHALF 3
#define QCCFILTER_SYMMETRIC_EXTENSION 0
#define QCCFILTER_PERIODIC_EXTENSION 1
#define QCCFILTER_SAMESAMPLING 0
#define QCCFILTER_SUBSAMPLEEVEN 1
#define QCCFILTER_SUBSAMPLEODD 2
#define QCCFILTER_UPSAMPLEEVEN 3
#define QCCFILTER_UPSAMPLEODD 4

typedef struct
{
  int causality;
  int length;
  QccVector coefficients;
} QccFilter;

int QccFilterInitialize(QccFilter *filter)
{
  if (filter == NULL)
    return(0);
  
  filter->causality = QCCFILTER_CAUSAL;
  filter->length = 0;
  filter->coefficients = NULL;

  return(0);
}


int QccFilterAlloc(QccFilter *filter)
{
  if (filter == NULL)
    return(0);
  
  if ((filter->coefficients == NULL) &&
      (filter->length > 0))
    if ((filter->coefficients = QccVectorAlloc(filter->length)) == NULL)
      {
        printf("(QccFilterAlloc): Error calling QccVectorAlloc()");
        return(1);
      }
  
  return(0);
}

static int QccFilterVectorPeriodicExtension(const QccVector input_signal,
                                            QccVector output_signal,
                                            int length,
                                            const QccFilter *filter)
{
  int input_index1;
  int input_index2;
  int output_index;
  int filter_index;


  switch (filter->causality)
    {

    case QCCFILTER_SYMMETRICHALF:
      for (output_index = 0; output_index < length; output_index++)
        for (filter_index = 0; filter_index < filter->length; filter_index++)
          {
            input_index1 = 
              QccMathModulus(output_index - filter_index,
                             length);
            input_index2 = 
              QccMathModulus(output_index + filter_index + 1,
                             length);
            
            output_signal[output_index] +=
              (input_signal[input_index1] + input_signal[input_index2]) * 
              filter->coefficients[filter_index];
          }
      break;
      
    default:
      printf("(QccFilterVectorPeriodicExtension): Undefined filter causality (%d)",
                         filter->causality);
      return(1);
    }
  
  return(0);
}


int QccFilterCalcSymmetricExtension(int index, int length, int symmetry)
{
  if (symmetry == QCCFILTER_SYMMETRICWHOLE)
    {
      if (length > 2)
        {
          if ((index < 0) || (index >= (2*length - 2)))
            index = QccMathModulus(index, 2*length - 2);
          if (index >= length)
            index = 2*length - index - 2;
        }
      else
        if ((index < 0) || (index >= length))
          index = QccMathModulus(index, length);
    }
  else
    {
      if (length > 1)
        {
          if ((index < 0) || (index >= 2*length))
            index = QccMathModulus(index, 2*length);
          if (index >= length)
            index = 2*length - index - 1;
        }
      else
        index = 0;
    }
  return(index);
}


static int QccFilterVectorSymmetricExtension(const QccVector input_signal,
                                             QccVector output_signal,
                                             int length,
                                             const QccFilter *filter)
{
  int input_index1;
  int input_index2;
  int output_index;
  int filter_index;

  switch (filter->causality)
    {
    
    case QCCFILTER_SYMMETRICHALF:
      for (output_index = 0; output_index < length; output_index++)
        for (filter_index = 0; filter_index < filter->length; filter_index++)
          {
            input_index1 = 
              QccFilterCalcSymmetricExtension(output_index - filter_index,
                                              length,
                                              QCCFILTER_SYMMETRICHALF);
            input_index2 = 
              QccFilterCalcSymmetricExtension(output_index + filter_index + 1,
                                              length,
                                              QCCFILTER_SYMMETRICHALF);
            
            output_signal[output_index] +=
              (input_signal[input_index1] + input_signal[input_index2]) * 
              filter->coefficients[filter_index];
        }
      break;
      
    default:
      printf("(QccFilterVectorSymmetricExtension): Undefined filter causality (%d)",
                         filter->causality);
      return(1);
    }
  
  return(0);
}


int QccFilterVector(const QccVector input_signal,
                    QccVector output_signal,
                    int length,
                    const QccFilter *filter,
                    int boundary_extension)
{
  if (input_signal == NULL)
    return(0);
  if (output_signal == NULL)
    return(0);
  if (filter == NULL)
    return(0);
  if (!(filter->length) || (filter->coefficients == NULL))
    return(0);
  
  switch (boundary_extension)
    {
    case QCCFILTER_SYMMETRIC_EXTENSION:
      if (QccFilterVectorSymmetricExtension(input_signal,
                                            output_signal,
                                            length,
                                            filter))
        {
          printf("(QccFilterVector): Error calling QccFilterVectorPeriodicExtension()");
          return(1);
        }
      break;
      
    case QCCFILTER_PERIODIC_EXTENSION:
      if (QccFilterVectorPeriodicExtension(input_signal,
                                           output_signal,
                                           length,
                                           filter))
        {
          printf("(QccFilterVector): Error calling QccFilterVectorSymmetricExtension()");
          return(1);
        }
      break;
      
    default:
      printf("(QccFilterVector): Undefined boundary extension (%d)",
                         boundary_extension);
      return(1);
    }
  
  return(0);
}

void QccFilterFree(QccFilter *filter)
{
  if (filter == NULL)
    return;
  
  QccVectorFree(filter->coefficients);
  filter->coefficients = NULL;
}



int QccFilterRead(FILE *infile,
                  QccFilter *filter)
{
  int index;
  
  if (infile == NULL)
    return(0);
  if (filter == NULL)
    return(0);
  
  fscanf(infile, "%d", &(filter->causality));
  if (ferror(infile) || feof(infile))
    goto Error;
  
  if (QccFileSkipWhiteSpace(infile, 0))
    goto Error;
  
  fscanf(infile, "%d", &(filter->length));
  if (ferror(infile) || feof(infile))
    goto Error;
  
  if (QccFileSkipWhiteSpace(infile, 0))
    goto Error;
  
  if (QccFilterAlloc(filter))
    {
      printf("(QccFilterRead): Error calling QccFilterAlloc()");
      return(1);
    }
  
  for (index = 0; index < filter->length; index++)
    {
      fscanf(infile, "%lf",
             &(filter->coefficients[index]));
      if (ferror(infile) || feof(infile))
        {
          QccFilterFree(filter);
          goto Error;
        }
    }
  
  return(0);
  
 Error:
  printf("(QccFilterRead): Error reading filter");
  return(1);
  
}


/*  huffman_codeword.c  */
#define QCCENTHUFFMAN_MAXCODEWORDLEN 31
typedef struct
{
  int length;
  int codeword;
} QccENTHuffmanCodeword;

/*  huffman_table.c  */
#define QCCENTHUFFMAN_DECODETABLE 0
#define QCCENTHUFFMAN_ENCODETABLE 1
#define QCCENTHUFFMAN_MAXSYMBOL 100000
#define QCCENTHUFFMANTABLE_MAGICNUM "HUF"
typedef struct
{
  int symbol;
  QccENTHuffmanCodeword codeword;
} QccENTHuffmanTableEntry;
typedef struct
{
  QccString filename;
  QccString magic_num;
  int major_version;
  int minor_version;
  int table_type;
  int table_length;
  QccENTHuffmanTableEntry *table;
  int *num_codewords_list;
  int num_codewords_list_length;
  int *symbol_list;
  int symbol_list_length;
  int codeword_max[QCCENTHUFFMAN_MAXCODEWORDLEN];
  int codeword_min[QCCENTHUFFMAN_MAXCODEWORDLEN];
  int codeword_ptr[QCCENTHUFFMAN_MAXCODEWORDLEN];
} QccENTHuffmanTable;


/*image_component.c*/

#define QCCIMGIMAGECOMPONENT_MAGICNUM "ICP"
typedef QccMatrix QccIMGImageArray;
typedef struct
{
  QccString filename;
  QccString magic_num;
  int major_version;
  int minor_version;
  int num_rows;
  int num_cols;
  double min_val;
  double max_val;
  QccIMGImageArray image;
} QccIMGImageComponent;



static QccIMGImageArray QccIMGImageArrayAlloc(int num_rows, int num_cols)
{
  QccIMGImageArray new_array;
  int row;
  
  if ((num_rows <= 0) && (num_cols <= 0))
    return(NULL);
  
  if ((new_array = (double **)malloc(sizeof(double *)*num_rows)) == NULL)
    {
      printf("(QccIMGImageArrayAlloc): Error allocating memory");
      return(NULL);
    }
  
  for (row = 0; row < num_rows; row++)
    if ((new_array[row] = (double *)malloc(sizeof(double)*num_cols)) == NULL)
      {
        printf("(QccIMGImageArrayAlloc): Error allocating memory");
        return(NULL);
      }
  
  return(new_array);
}


static void QccIMGImageArrayFree(QccIMGImageArray image_array, int num_rows)
{
  int row;
  
  if (image_array != NULL)
    {
      for (row = 0; row < num_rows; row++)
        if (image_array[row] != NULL)
          QccFree(image_array[row]);
      QccFree(image_array);
    }
}


int QccIMGImageComponentInitialize(QccIMGImageComponent *image_component)
{
  if (image_component == NULL)
    return(0);

  QccStringMakeNull(image_component->filename);
  QccStringCopy(image_component->magic_num, QCCIMGIMAGECOMPONENT_MAGICNUM);

  image_component->num_rows = 0;
  image_component->num_cols = 0;
  image_component->min_val = 0.0;
  image_component->max_val = 0.0;
  image_component->image = NULL;

  return(0);
}



int QccIMGImageComponentAlloc(QccIMGImageComponent *image_component)
{
  if (image_component == NULL)
    return(0);
  
  if (image_component->image == NULL)
    {
      if ((image_component->num_rows > 0) ||
          (image_component->num_cols > 0))
        {
          if ((image_component->image = 
               QccIMGImageArrayAlloc(image_component->num_rows, 
                                     image_component->num_cols)) == NULL)
            {
              printf("(QccIMGImageComponentAlloc): Error calling QccIMGImageArrayAlloc()");
              return(1);
            }
        }
      else
        image_component->image = NULL;
    }
  
  return(0);
}


void QccIMGImageComponentFree(QccIMGImageComponent *image_component)
{
  if (image_component == NULL)
    return;
  
  if (image_component->image != NULL)
    {
      QccIMGImageArrayFree(image_component->image, 
                           image_component->num_rows);
      image_component->image = NULL;
    }
}





int QccIMGImageComponentSetMin(QccIMGImageComponent *image_component)
{
  double min = MAXDOUBLE;
  int row, col;
  
  if (image_component == NULL)
    return(0);
  if (image_component->image == NULL)
    return(0);
  
  for (row = 0; row < image_component->num_rows; row++)
    for (col = 0; col < image_component->num_cols; col++)
      if (image_component->image[row][col] < min)
        min = image_component->image[row][col];
  
  image_component->min_val = min;
  
  return(0);
}


int QccIMGImageComponentSetMax(QccIMGImageComponent *image_component)
{
  double max = -MAXDOUBLE;
  int row, col;
  
  if (image_component == NULL)
    return(0);
  if (image_component->image == NULL)
    return(0);
  
  for (row = 0; row < image_component->num_rows; row++)
    for (col = 0; col < image_component->num_cols; col++)
      if (image_component->image[row][col] > max)
        max = image_component->image[row][col];
  
  image_component->max_val = max;
  
  return(0);
}

int QccIMGImageComponentSetMaxMin(QccIMGImageComponent *image_component)
{
  if (QccIMGImageComponentSetMax(image_component))
    {
      printf("(QccIMGImageComponentSetMaxMin): Error calling QccIMGImageComponentSetMax()");
      return(1);
    }

  if (QccIMGImageComponentSetMin(image_component))
    {
      printf("(QccIMGImageComponentSetMaxMin): Error calling QccIMGImageComponentSetMin()");
      return(1);
    }

  return(0);
}

int QccIMGImageComponentSubtract(const QccIMGImageComponent *image_component1,
                                 const QccIMGImageComponent *image_component2,
                                 QccIMGImageComponent *image_component3)
{
  int row, col;

  if (image_component1 == NULL)
    return(0);
  if (image_component2 == NULL)
    return(0);
  if (image_component3 == NULL)
    return(0);

  if (image_component1->image == NULL)
    return(0);
  if (image_component2->image == NULL)
    return(0);
  if (image_component3->image == NULL)
    return(0);

  if ((image_component1->num_rows != image_component2->num_rows) ||
      (image_component1->num_cols != image_component2->num_cols) ||
      (image_component1->num_rows != image_component3->num_rows) ||
      (image_component1->num_cols != image_component3->num_cols))
    {
      printf("(QccIMGImageComponentSubtract): Image components must have same number of rows and columns");
      return(1);
    }

  for (row = 0; row < image_component1->num_rows; row++)
    for (col = 0; col < image_component1->num_cols; col++)
      image_component3->image[row][col] =
        image_component1->image[row][col] -
        image_component2->image[row][col];

  if (QccIMGImageComponentSetMaxMin(image_component3))
    {
      printf("(QccIMGImageComponentSubtract): Error calling QccIMGImageComponentSetMaxMin()");
      return(1);
    }

  return(0);
}


int QccIMGImageComponentAbsoluteValue(QccIMGImageComponent *image_component)
{
  int row, col;

  if (image_component == NULL)
    return(0);
  if (image_component->image == NULL)
    return(0);

  for (row = 0; row < image_component->num_rows; row++)
    for (col = 0; col < image_component->num_cols; col++)
      image_component->image[row][col] =
        fabs(image_component->image[row][col]);

  if (QccIMGImageComponentSetMaxMin(image_component))
    {
      printf("(QccIMGImageComponentAbsoluteValue): Error calling QccIMGImageComponentSetMaxMin()");
      return(1);
    }

  return(0);
}



static int QccIMGImageComponentReadHeader(FILE *infile, 
                                          QccIMGImageComponent
                                          *image_component)
{
  if ((infile == NULL) || (image_component == NULL))
    return(0);
  
  if (QccFileReadMagicNumber(infile,
                             image_component->magic_num,
                             &image_component->major_version,
                             &image_component->minor_version))
    {
      printf("(QccIMGImageComponentReadHeader): Error reading magic number in %s",
                         image_component->filename);
      return(1);
    }
  
  if (strcmp(image_component->magic_num, QCCIMGIMAGECOMPONENT_MAGICNUM))
    {
      printf("(QccIMGImageComponentReadHeader): %s is not of image_component (%s) type",
                         image_component->filename,
                         QCCIMGIMAGECOMPONENT_MAGICNUM);
      return(1);
    }
  
  fscanf(infile, "%d", &(image_component->num_cols));
  if (ferror(infile) || feof(infile))
    {
      printf("(QccIMGImageComponentReadHeader): Error reading number of columns in %s",
                         image_component->filename);
      return(1);
    }
  
  if (QccFileSkipWhiteSpace(infile, 0))
    {
      printf("(QccIMGImageComponentReadHeader): Error reading number of columns in %s",
                         image_component->filename);
      return(1);
    }
  
  fscanf(infile, "%d", &(image_component->num_rows));
  if (ferror(infile) || feof(infile))
    {
      printf("(QccIMGImageComponentReadHeader): Error reading number of rows in %s",
                         image_component->filename);
      return(1);
    }
  
  if (QccFileSkipWhiteSpace(infile, 0))
    {
      printf("(QccIMGImageComponentReadHeader): Error reading number of rows in %s",
                         image_component->filename);
      return(1);
    }
  
  fscanf(infile, "%lf", &(image_component->min_val));
  if (ferror(infile) || feof(infile))
    {
      printf("(QccIMGImageComponentReadHeader): Error reading number of rows in %s",
                         image_component->filename);
      return(1);
    }
  
  if (QccFileSkipWhiteSpace(infile, 0))
    {
      printf("(QccIMGImageComponentReadHeader): Error reading number of rows in %s",
                         image_component->filename);
      return(1);
    }
  
  fscanf(infile, "%lf%*1[\n]", &(image_component->max_val));
  if (ferror(infile) || feof(infile))
    {
      printf("(QccIMGImageComponentReadHeader): Error reading number of rows in %s",
                         image_component->filename);
      return(1);
    }
  
  return(0);
}


static int QccIMGImageComponentReadData(FILE *infile, 
                                        QccIMGImageComponent *image_component)
{
  int row, col;
  
  if ((infile == NULL) ||
      (image_component == NULL))
    return(0);
  
  for (row = 0; row < image_component->num_rows; row++)
    for (col = 0; col < image_component->num_cols; col++)
      if (QccFileReadDouble(infile,
                            &(image_component->image[row][col])))
        {
          printf("(QccIMGImageComponentReadData): Error calling QccFileReadDouble()",
                             image_component->filename);
          return(1);
        }
  
  return(0);
}


int QccIMGImageComponentRead(QccIMGImageComponent *image_component)
{
  FILE *infile = NULL;
  
  if (image_component == NULL)
    return(0);
  
  if ((infile = QccFileOpen(image_component->filename, "r")) == NULL)
    {
      printf("(QccIMGImageComponentRead): Error opening %s for reading",
                         image_component->filename);
      return(1);
    }
  
  if (QccIMGImageComponentReadHeader(infile, image_component))
    {
      printf("(QccIMGImageComponentRead): Error calling QccIMGImageComponentReadHeader()");
      return(1);
    }
  
  if (QccIMGImageComponentAlloc(image_component))
    {
      printf("(QccIMGImageComponentRead): Error calling QccIMGImageComponentAlloc()");
      return(1);
    }
  
  if (QccIMGImageComponentReadData(infile, image_component))
    {
      printf("(QccIMGImageComponentRead): Error calling QccIMGImageComponentReadData()");
      return(1);
    }
  
  QccFileClose(infile);
  return(0);
}


static int QccIMGImageComponentWriteHeader(FILE *outfile, 
                                           const QccIMGImageComponent 
                                           *image_component)
{
  if ((outfile == NULL) || (image_component == NULL))
    return(0);
  
  if (QccFileWriteMagicNumber(outfile, QCCIMGIMAGECOMPONENT_MAGICNUM))
    goto Error;
  
  fprintf(outfile, "%d %d\n",
          image_component->num_cols,
          image_component->num_rows);
  if (ferror(outfile))
    goto Error;
  
  fprintf(outfile, "% 16.9e % 16.9e\n",
          image_component->min_val,
          image_component->max_val);
  if (ferror(outfile))
    goto Error;
  
  return(0);
  
 Error:
  printf("(QccIMGImageComponentWriteHeader): Error writing header to %s",
                     image_component->filename);
  return(1);
}


static int QccIMGImageComponentWriteData(FILE *outfile,
                                         const QccIMGImageComponent
                                         *image_component)
{
  int row, col;
  
  if ((image_component == NULL) ||
      (outfile == NULL))
    return(0);
  
  for (row = 0; row < image_component->num_rows; row++)
    for (col = 0; col < image_component->num_cols; col++)
      if (QccFileWriteDouble(outfile,
                             image_component->image[row][col]))
        {
          printf("(QccIMGImageComponentWriteData): Error calling QccFileWriteDouble()");
          return(1);
        }
  
  return(0);
}


int QccIMGImageComponentWrite(const QccIMGImageComponent *image_component)
{
  FILE *outfile;
  
  if (image_component == NULL)
    return(0);
  
  if ((outfile = QccFileOpen(image_component->filename, "w")) == NULL)
    {
      printf("(QccIMGImageComponentWrite): Error opening %s for writing",
                         image_component->filename);
      return(1);
    }
  
  if (QccIMGImageComponentWriteHeader(outfile, image_component))
    {
      printf("(QccIMGImageComponentWrite): Error calling QccIMGImageComponentWriteHeader()");
      return(1);
    }
  if (QccIMGImageComponentWriteData(outfile, image_component))
    {
      printf("(QccIMGImageComponentWrite): Error calling QccIMGImageComponentWriteData()");
      return(1);
    }
  
  QccFileClose(outfile);
  return(0);
}

int QccIMGImageComponentCopy(QccIMGImageComponent *image_component1,
                             const QccIMGImageComponent *image_component2)
{
  int row, col;
  
  if ((image_component1 == NULL) ||
      (image_component2 == NULL))
    return(0);
  
  if ((image_component2->image == NULL) ||
      (image_component2->num_rows <= 0) ||
      (image_component2->num_cols <= 0))
    return(0);

  if (image_component1->image == NULL)
    {
      image_component1->num_rows = image_component2->num_rows;
      image_component1->num_cols = image_component2->num_cols;
      if (QccIMGImageComponentAlloc(image_component1))
        {
          printf("(QccIMGImageComponentCopy): Error calling QccIMGImageComponentAlloc()");
          return(1);
        }
    }
  else
    {
      if ((image_component1->num_rows != image_component2->num_rows) ||
          (image_component1->num_cols != image_component2->num_cols))
        {
          printf("(QccIMGImageComponentCopy): Image-component arrays have different sizes");
          return(1);
        }
    }

  for (row = 0; row < image_component1->num_rows; row++)
    for (col = 0; col < image_component1->num_cols; col++)
      image_component1->image[row][col] = 
        image_component2->image[row][col];
  
  if (QccIMGImageComponentSetMaxMin(image_component1))
    {
      printf("(QccIMGImageComponentCopy): Error calling QccIMGImageComponentSetMaxMin()");
      return(1);
    }

  return(0);
  
}




int QccIMGImageComponentInterpolateBilinear(const QccIMGImageComponent
                                            *image_component1,
                                            QccIMGImageComponent
                                            *image_component2)
{
  int row1, col1;
  int row2, col2;

  if (image_component1 == NULL)
    return(0);
  if (image_component2 == NULL)
    return(0);

  if (image_component1->image == NULL)
    return(0);
  if (image_component2->image == NULL)
    return(0);

  if ((image_component2->num_rows != 2 * image_component1->num_rows) ||
      (image_component2->num_cols != 2 * image_component1->num_cols))
    {
      printf("(QccIMGImageComponentInterpolateBilinear): Destination image size must be twice that of source image");
      return(1);
    }
      
  for (row1 = 0, row2 = 0; row1 < image_component1->num_rows;
       row1++, row2 += 2)
    for (col1 = 0, col2 = 0; col1 < image_component1->num_cols;
         col1++, col2 += 2)
      image_component2->image[row2][col2] =
        image_component1->image[row1][col1];

  for (row2 = 1; row2 < image_component2->num_rows - 1; row2 += 2)
    for (col2 = 0; col2 < image_component2->num_cols; col2 += 2)
      image_component2->image[row2][col2] =
        (image_component2->image[row2 - 1][col2] +
         image_component2->image[row2 + 1][col2]) / 2;

  for (col2 = 1; col2 < image_component2->num_cols - 1; col2 +=2)
    for (row2 = 0; row2 < image_component2->num_rows; row2++)
      image_component2->image[row2][col2] =
        (image_component2->image[row2][col2 - 1] +
         image_component2->image[row2][col2 + 1]) / 2;

  for (col2 = 0; col2 < image_component2->num_cols; col2++)
    image_component2->image[image_component2->num_rows - 1][col2] =
      image_component2->image[image_component2->num_rows - 2][col2];
  for (row2 = 0; row2 < image_component2->num_rows; row2++)
    image_component2->image[row2][image_component2->num_cols - 1] = 
      image_component2->image[row2][image_component2->num_cols - 2];

  if (QccIMGImageComponentSetMaxMin(image_component2))
    {
      printf("(QccIMGImageComponentInterpolateBilinear): Error calling QccIMGImageComponentSetMaxMin()");
      return(1);
    }

  return(0);
}



int QccIMGImageComponentInterpolateFilter(const QccIMGImageComponent
                                          *image_component1,
                                          QccIMGImageComponent *image_component2,
                                          const QccFilter *filter)
{
  int return_value;
  int row1, col1;
  int row2, col2;
  QccVector input_vector = NULL;
  QccVector output_vector = NULL;

  if (image_component1 == NULL)
    return(0);
  if (image_component2 == NULL)
    return(0);
  if (filter == NULL)
    return(0);

  if (image_component1->image == NULL)
    return(0);
  if (image_component2->image == NULL)
    return(0);

  if (filter->causality != QCCFILTER_SYMMETRICHALF)
    {
      printf("(QccIMGImageComponentInterpolateFilter): Filter must be half-sample symmetric");
      goto Error;
    }

  if ((image_component2->num_rows != 2 * image_component1->num_rows) ||
      (image_component2->num_cols != 2 * image_component1->num_cols))
    {
      printf("(QccIMGImageComponentInterpolateFilter): Destination image size must be twice that of source image");
      goto Error;
    }
      
  if ((input_vector =
       QccVectorAlloc(QccMathMax(image_component2->num_rows,
                                 image_component2->num_cols))) == NULL)
    {
      printf("(QccIMGImageComponentInterpolateFilter): Error calling QccVectorAlloc()");
      goto Error;
    }
  if ((output_vector =
       QccVectorAlloc(QccMathMax(image_component2->num_rows,
                                 image_component2->num_cols))) == NULL)
    {
      printf("(QccIMGImageComponentInterpolateFilter): Error calling QccVectorAlloc()");
      goto Error;
    }
  
  for (row1 = 0, row2 = 0; row1 < image_component1->num_rows;
       row1++, row2 += 2)
    for (col1 = 0, col2 = 0; col1 < image_component1->num_cols;
         col1++, col2 += 2)
      image_component2->image[row2][col2] =
        image_component1->image[row1][col1];
  
  for (col1 = 0; col1 < image_component1->num_cols; col1++)
    {
      for (row1 = 0; row1 < image_component1->num_rows; row1++)
        input_vector[row1] = image_component1->image[row1][col1];

      QccVectorZero(output_vector, image_component1->num_rows);

      if (QccFilterVector(input_vector,
                          output_vector,
                          image_component1->num_rows,
                          filter,
                          QCCFILTER_SYMMETRIC_EXTENSION))
        {
          printf("(QccIMGImageComponentInterpolateFilter): Error calling QccFilterVector()");
          goto Error;
        }

      for (row2 = 1; row2 < image_component2->num_rows - 1; row2 += 2)
        image_component2->image[row2][col1 * 2] =
          output_vector[row2 / 2];

      image_component2->image[image_component2->num_rows - 1][col1 * 2] =
        image_component2->image[image_component2->num_rows - 2][col1 * 2];
    }
  
  for (row2 = 0; row2 < image_component2->num_rows; row2++)
    {
      for (col1 = 0; col1 < image_component1->num_cols; col1++)
        input_vector[col1] =
          image_component2->image[row2][col1 * 2];

      QccVectorZero(output_vector, image_component1->num_cols);

      if (QccFilterVector(input_vector,
                          output_vector,
                          image_component1->num_cols,
                          filter,
                          QCCFILTER_SYMMETRIC_EXTENSION))
        {
          printf("(QccIMGImageComponentInterpolateFilter): Error calling QccFilterVector()");
          goto Error;
        }
      
      for (col2 = 1; col2 < image_component2->num_cols - 1; col2 +=2)
        image_component2->image[row2][col2] =
          output_vector[col2 / 2];

      image_component2->image[row2][image_component2->num_cols - 1] =
        image_component2->image[row2][image_component2->num_cols - 2];
    }
  
  if (QccIMGImageComponentSetMaxMin(image_component2))
    {
      printf("(QccIMGImageComponentInterpolateFilter): Error calling QccIMGImageComponentSetMaxMin()");
      goto Error;
    }
  
  return_value = 0;
  goto Return;
 Error:
  return_value = 1;
 Return:
  QccVectorFree(input_vector);
  QccVectorFree(output_vector);
  return(return_value);
}


/*  parse.c  */

#define QCCPARSEMAXNUMARGUMENTS 1024
#define QCCPARSESTRINGLEN QCCSTRINGLEN


static void QccParseCatChar(char *str, char ch)
{
  int length;

  length = strlen(str);

  if (length > (QCCPARSESTRINGLEN - 1))
    return;

  str[length] = ch;
  str[length + 1] = '\0';
}


static const char *QccParseGetFormat(char *fmt, const char *pnt,
                                     int *multiple_arguments)
{
  int i = 0;

  *multiple_arguments = 0;
  while ((*pnt != '\0') && (*pnt != ':') && (*pnt != ']') && (*pnt != '[') && 
         (*pnt != ' ') && (*pnt != '\n') && (*pnt != '\t') && (*pnt != '\r'))
    {
      if (*pnt == '*')
        *multiple_arguments = 1;
      else
        fmt[i++] = *pnt;
      pnt = &(pnt[1]);
    }
  fmt[i] = '\0';

  return(pnt);
}


static char QccParseGetType(const char *pnt, int *multiple_arguments)
{
  char type = 'd';

  *multiple_arguments = 0;

  while ((*pnt != '\0') && (*pnt != ':') && (*pnt != ']') && (*pnt != '[') && 
         (*pnt != ' ') && (*pnt != '\n') && (*pnt != '\t') && (*pnt != '\r'))
    {
      if (pnt[0] == '*')
        *multiple_arguments = 1;
      else
        type = pnt[0];
      pnt = &(pnt[1]);
    }

  switch (type)
    {
      /*  Unsigned integer  */
    case 'u':
      type = 'u';
      break;

      /*  Floating point  */
    case 'e':
    case 'g':
    case 'f':
      type = 'f';
      break;

      /*  Character string  */
    case 's':
      type = 's';
      break;

      /*  Integer  */
    default:
      type = 'd';
      break;
    }

  if (*multiple_arguments)
    {
      if (strchr(pnt, ']') != NULL)
        {
          printf("(QccParseParameters): Multiple argument designation can only be used for non-optional parameters");
          return(0);
        }
      if (strchr(pnt, ' ') != NULL)
        {
          printf("(QccParseParameters): Multiple argument designation must be last argument");
          return(0);
        }
    }

  return (type);
}


static char *QccParseFindPosition(const char *format, char *arg, int *pos)
{
  const char *format_orig;
  char *pnt;
  const char *tmp;
  char sw[QCCPARSESTRINGLEN + 1];
  int i, done = 0;;

  format_orig = format;

  do
    {
      sw[0] = '\0';
      /*  Find location of switch  */
      pnt = strstr(format, arg);
      if (pnt != NULL)
        {
          i = 0;
          /*  Extract full switch from format prototype  */
          while ((pnt[i] != '\0') && (pnt[i] != ':') && 
                 (pnt[i] != ']') && (pnt[i] != '[') && (pnt[i] != ' ') && 
                 (pnt[i] != '\n') && (pnt[i] != '\t') && (pnt[i] != '\r'))
            QccParseCatChar(sw, pnt[i++]);
          if (pnt[i] == '\0')
            pnt = NULL;
          /*  Make sure switches match exactly  */
          if (!strcmp(sw, arg))
            done = 1;
          else
            format = &pnt[i];
        }
    }
  while ((pnt != NULL) && (!done));

  /*  Count number of pointer references in prototype to current pos  */
  *pos = 0;
  tmp = format_orig;
  if (pnt != NULL)
    while (tmp != pnt)
      {
        if (tmp[0] == '%')
          {
            (*pos)++;
          }
        tmp = &tmp[1];
      }
  return(pnt);
}


static int QccParseReadParameter(int narg,
                                 char *arg[], 
                                 char *fmt,
                                 int pos,
                                 int *cindex,
                                 int multiple_arguments,
                                 void **parse_pointer,
                                 char *parse_type)
{
  int val = 0;
  int *num_args;
  int **d_ptr;
  unsigned int **u_ptr;
  float **f_ptr;
  QccString **s_ptr;

  if (!multiple_arguments)
    /*  Use type to cast the pointer appropriately  */
    switch (parse_type[pos])
      {
      case 'f':
        val = sscanf(arg[*cindex], fmt, (float *) parse_pointer[pos]);
        break;
      case 'u':
        val = sscanf(arg[*cindex], fmt, (unsigned int *) parse_pointer[pos]);
        break;
      case 's':
        QccConvertToQccString((char *)(parse_pointer[pos]),
                              arg[*cindex]);
        val = strlen((char *) parse_pointer[pos]);
        break;
      default:
        val = sscanf(arg[*cindex], fmt, (int *) parse_pointer[pos]);
        break;
      }
  else
    {
      num_args = (int *)parse_pointer[pos];
      *num_args = 1;
      for ( ; *cindex < narg; (*cindex)++, 
              *num_args += 1)
        switch (parse_type[pos])
          {
          case 'f':
            f_ptr = (float **)parse_pointer[pos + 1];
            if (*num_args == 1)
              {
                if ((*f_ptr = (float *)malloc(sizeof(float))) == NULL)
                  return(0);
              }
            else
              if ((*f_ptr = 
                   (float *)realloc(*f_ptr, sizeof(float)*(*num_args))) == 
                  NULL)
                return(0);
            val = sscanf(arg[*cindex], fmt, &(*f_ptr)[*num_args - 1]);
            if (val != 1)
              return(val);
            break;
          case 'u':
            u_ptr = (unsigned int **)parse_pointer[pos + 1];
            if (*num_args == 1)
              {
                if ((*u_ptr = 
                     (unsigned int *)malloc(sizeof(unsigned int))) == NULL)
                  return(0);
              }
            else
              if ((*u_ptr =
                   (unsigned int *)realloc(*u_ptr, 
                                           sizeof(unsigned int)*(*num_args))) 
                  == NULL)
                return(0);
            val = sscanf(arg[*cindex], fmt, &(*u_ptr)[*num_args - 1]);
            if (val != 1)
              return(val);
            break;
          case 's':
            s_ptr = (QccString **)parse_pointer[pos + 1];
            if (*num_args == 1)
              {
                if ((*s_ptr = (QccString *)malloc(sizeof(QccString))) == NULL)
                  return(0);
              }
            else
              if ((*s_ptr = 
                   (QccString *)realloc(*(s_ptr),
                                        sizeof(QccString) * (*num_args))) ==
                  NULL)
                return(0);
            QccConvertToQccString((*s_ptr)[*num_args - 1], arg[*cindex]);
            val = strlen(arg[*cindex]);
            if (val < 1)
              return(val);
            break;
          case 'd':
            d_ptr = (int **)parse_pointer[pos + 1];
            if (*num_args == 1)
              {
                if ((*d_ptr = (int *)malloc(sizeof(int))) == NULL)
                  return(0);
              }
            else
              if ((*d_ptr =
                   (int *)realloc(*d_ptr, sizeof(int)*(*num_args))) == NULL)
                return(0);
            val = sscanf(arg[*cindex], fmt, &((*d_ptr)[*num_args - 1]));
            if (val != 1)
              return(val);
            break;
          }
      *num_args -= 1;
    }

  return(val);
}


static const char *QccParseStartOfMandatoryParameters(const char *format,
                                                      int *pos)
{
  int i, done;
  const char *format2, *tmp;

  format2 = tmp = format;
  *pos = 0;

  /*  Find first ']' from end of string  */
  /*  This is last of optional parameters  */
  i = strlen(format2) - 1;
  done = 0;
  while ((i > 0) && (!done))
    {
      done = (format2[i] == ']');
      if (!done)
        i--;
    }

  /*  Scan until we find '%'  */
  /*  This is first mandatory parameter  */
  format2 = &format2[i];
  while((format2[0] != '\0') && (format2[0] != '%'))
    format2 = &format2[1];

  /*  Get count of the pointer reference  */
  while (tmp != format2)
    {
      if (tmp[0] == '%')
        (*pos)++;
      tmp = &tmp[1];
    }

  return(format2);
}


static void QccParsePrintUsage(const char *format)
{
  QccString program_name;
  QccString usage_string;

  QccStringMakeNull(program_name);
  QccStringMakeNull(usage_string);

  QccGetProgramName(program_name);

  sprintf(usage_string, "Usage: %s ", program_name);

  while (format[0] !=  '\0')
    {
      /*  Skip to label  */
      if (format[0] == '%')
        {
          while ((format[0] != ':') && (format[0] != ']') && 
                 (format[0] != '[') && (format[0] != ' ') && 
                 (format[0] != '\n') && (format[0] != '\t') && 
                 (format[0] != '\r'))
            format = &format[1];

          if (format[0] == ':')
            format = &format[1];
        }
      /*  Print labels and switches  */
      if (format[0] != '\0')
        {
          QccParseCatChar(usage_string, format[0]);
          format = &format[1];
        }
    }

  printf(usage_string);
}


int QccParseParametersVA(int argc, char *argv[], const char *format,
                         va_list ap)
{
  int done, switch_done;
  const char *format2;
  const char *pnt;
  char fmt[QCCPARSESTRINGLEN + 1];
  int cindex, i, pos;
  int multiple_arguments;
  int return_value = 0;
  void **parse_pointer = NULL;
  char *parse_type = NULL;

  if ((parse_pointer = 
       (void **)calloc(QCCPARSEMAXNUMARGUMENTS, sizeof(void *))) == NULL)
    {    
      printf("(QccParseParameters): Error allocating memory");
      goto Error;
    }

  if ((parse_type = 
       (char *)calloc(QCCPARSEMAXNUMARGUMENTS, sizeof(char))) == NULL)
    {    
      printf("(QccParseParameters): Error allocating memory");
      goto Error;
    }

  /*  Get pointers  */
  i = 0;
  pnt = &format[0];
  while ((i < QCCPARSEMAXNUMARGUMENTS) && (pnt[0] != '\0'))
    {
      if (pnt[0] == '%')
        {
          parse_type[i] = QccParseGetType(pnt, &multiple_arguments);
          if (multiple_arguments)
            switch (parse_type[i])
              {
              case 'f':
                parse_pointer[i++] = (void *) va_arg(ap, int *);
                parse_pointer[i] = (void *) va_arg(ap, float **);
                *((float **)parse_pointer[i++]) = NULL;
                break;
              case 'u':
                parse_pointer[i++] = (void *) va_arg(ap, int *);
                parse_pointer[i] = (void *) va_arg(ap, unsigned int **);
                *((unsigned int **)parse_pointer[i++]) = NULL;
                break;
              case 's':
                parse_pointer[i++] = (void *) va_arg(ap, int *);
                parse_pointer[i] = (void *) va_arg(ap, char ***);
                *((unsigned char ***)parse_pointer[i++]) = NULL;
                break;
              case 'd':
                parse_pointer[i++] = (void *) va_arg(ap, int *);
                parse_pointer[i] = (void *) va_arg(ap, int **);
                *((int **)parse_pointer[i++]) = NULL;
                break;
              }
          else
            switch (parse_type[i])
              {
              case 'f':
                parse_pointer[i++] = (void *) va_arg(ap, float *);
                break;
              case 'u':
                parse_pointer[i++] = (void *) va_arg(ap, unsigned int *);
                break;
              case 's':
                parse_pointer[i++] = (void *) va_arg(ap, char *);
                break;
              case 'd':
                parse_pointer[i++] = (void *) va_arg(ap, int *);
                break;
              default:
                goto QccParseErrorInUsageStringReturn;
              }
        }
      pnt = &pnt[1];
    }
  if (i >= QCCPARSEMAXNUMARGUMENTS)
    {
      printf("(QccParseParameters): Too many arguments");
      goto QccParseErrorInUsageStringReturn;
    }

  /*  Process all switches first  */
  cindex = 1;
  switch_done = 0;
  while ((cindex < argc) && !switch_done)
    {
      /*  If not a switch then done  */
      if (argv[cindex][0] != '-')
        switch_done = 1;
      else
        {
          /*  Find switch in prototype  */
          pnt = QccParseFindPosition(format, argv[cindex], &pos);
          if (pnt == NULL)
            goto QccParseParseErrorReturn;
          cindex++;
          done = 0;
          while(!done)
            {
              /*  Pointer reference  */
              if (pnt[0] == '%')
                {
                  /*  Pointer flag, no associated value  */
                  if (pnt[1] == ':')
                    {
                      /*  Set the flag  */
                      *(int *)parse_pointer[pos++] = 1;
                      pnt = &pnt[2];
                    }
                  else
                    {
                      if (cindex >= argc)
                        goto QccParseParseErrorReturn;
                      /*  Get format for argument  */
                      pnt = QccParseGetFormat(fmt, pnt, &multiple_arguments);
                      /*  Get value of argument  */
                      if (!QccParseReadParameter(argc, argv,
                                                 fmt, pos++, &cindex, 0,
                                                 parse_pointer,
                                                 parse_type))
                        goto QccParseParseErrorReturn;
                      cindex++;
                    }
                }
              else
                {
                  /*  Keep searching till end of switch's parameters  */
                  if ((pnt[0] == '\0')||(pnt[0] == ']')||(pnt[0] == '['))
                    done = 1;
                  else
                    pnt = &pnt[1];
                }
            }
        }
    }

  /*  Now do all mandatory parameters  */
  format2 = QccParseStartOfMandatoryParameters(format, &pos);

  /*  Loop till end of string or out of arguments  */
  while ((format2[0] != '\0') && (cindex < argc))
    {
      /*  Pointer reference  */
      if (format2[0] == '%')
        {
          /*  Get format for argument  */
          format2 = QccParseGetFormat(fmt, format2, &multiple_arguments);
          /*  Get value of argument  */
          if (!QccParseReadParameter(argc, argv,
                                     fmt, pos++, &cindex, multiple_arguments,
                                     parse_pointer, parse_type))
            goto QccParseParseErrorReturn;
        }
      
      cindex++;
      /*  Loop to next pointer reference  */
      while((format2[0] != '\0') && (format2[0] != '%'))
        format2 = &format2[1];
    }

  /*  Check to see if all mandatory paramters have been processed  */
  if ((format2[0] != '\0') || (cindex < argc))
    goto QccParseParseErrorReturn;

  /* Normal error-free return */
  return_value = 0;
  goto Return;
  
 QccParseParseErrorReturn:
  QccParsePrintUsage(format);
  goto Error;

 QccParseErrorInUsageStringReturn:
  printf("(QccParseParameters): Error in usage string");
  goto Error;

 Error:
  return_value = 1;

 Return:
  if (parse_pointer != NULL)
    QccFree(parse_pointer);
  if (parse_type != NULL)
    QccFree(parse_type);
  return(return_value);
}


int QccParseParameters(int argc, char *argv[], const char *format, ...)
{
  int return_value;
  va_list ap;
  
  va_start(ap, format);
  
  return_value = QccParseParametersVA(argc, argv, format, ap);

  va_end(ap);

  return(return_value);
}


/*  motion_vectors.c  */
typedef struct
{
  QccENTHuffmanTable encode_table;
  QccENTHuffmanTable decode_table;
} QccVIDMotionVectorsTable;


int QccVIDMotionVectorsWriteFile(const QccIMGImageComponent
                                 *motion_vectors_horizontal,
                                 const QccIMGImageComponent
                                 *motion_vectors_vertical,
                                 const QccString filename,
                                 int frame_num)
{
  QccString current_filename;
  FILE *current_file;
  int row, col;

  if (motion_vectors_horizontal == NULL)
    return(0);
  if (motion_vectors_vertical == NULL)
    return(0);
  if (filename == NULL)
    return(0);
  if (QccStringNull(filename))
    return(0);

  if (frame_num < 0)
    {
      printf("(QccVIDMotionVectorsWriteFile): Invalid frame number (%d)",
                         frame_num);
      return(1);
    }

  if ((motion_vectors_horizontal->num_rows !=
       motion_vectors_vertical->num_rows) ||
      (motion_vectors_horizontal->num_cols !=
       motion_vectors_vertical->num_cols))
    {
      printf("(QccVIDMotionVectorsWriteFile): Motion-vector fields must have same size");
      return(1);
    }

  QccStringSprintf(current_filename, filename, frame_num);

  if ((current_file = QccFileOpen(current_filename, "w")) == NULL)
    {
      printf("(QccVIDMotionVectorsWriteFile): Error calling QccFileOpen()");
      return(1);
    }

  for (row = 0; row < motion_vectors_horizontal->num_rows; row++)
    for (col = 0; col < motion_vectors_horizontal->num_cols; col++)
      fprintf(current_file, "% 11.4f % 11.4f\n",
              motion_vectors_vertical->image[row][col],
              motion_vectors_horizontal->image[row][col]);

  QccFileClose(current_file);

  return(0);
}




/*image.c*/
#define QCCIMGTYPE_UNKNOWN 0
#define QCCIMGTYPE_PBM     1
#define QCCIMGTYPE_PGM     2
#define QCCIMGTYPE_PPM     3
#define QCCIMGTYPE_ICP     4
typedef struct
{
  int image_type;
  QccString filename;
  QccIMGImageComponent Y;
  QccIMGImageComponent U;
  QccIMGImageComponent V;
} QccIMGImage;



int QccIMGImageInitialize(QccIMGImage *image)
{
  if (image == NULL)
    return(0);
  
  image->image_type = QCCIMGTYPE_UNKNOWN;
  QccStringMakeNull(image->filename);
  QccIMGImageComponentInitialize(&(image->Y));
  QccIMGImageComponentInitialize(&(image->U));
  QccIMGImageComponentInitialize(&(image->V));
  
  return(0);
}


int QccIMGImageGetSize(const QccIMGImage *image,
                       int *num_rows, int *num_cols)
{
  if (image == NULL)
    return(0);
  
  if (num_rows != NULL)
    *num_rows = image->Y.num_rows;
  if (num_cols != NULL)
    *num_cols = image->Y.num_cols;
  
  return(0);
}


int QccIMGImageGetSizeYUV(const QccIMGImage *image,
                          int *num_rows_Y, int *num_cols_Y,
                          int *num_rows_U, int *num_cols_U,
                          int *num_rows_V, int *num_cols_V)
{
  if (image == NULL)
    return(0);
  
  if (num_rows_Y != NULL)
    *num_rows_Y = image->Y.num_rows;
  if (num_cols_Y != NULL)
    *num_cols_Y = image->Y.num_cols;
  
  if (num_rows_U != NULL)
    *num_rows_U = image->U.num_rows;
  if (num_cols_U != NULL)
    *num_cols_U = image->U.num_cols;
  
  if (num_rows_V != NULL)
    *num_rows_V = image->V.num_rows;
  if (num_cols_V != NULL)
    *num_cols_V = image->V.num_cols;
  
  return(0);
}


int QccIMGImageSetSize(QccIMGImage *image,
                       int num_rows, int num_cols)
{
  if (image == NULL)
    return(0);
  
  if (QccIMGImageColor(image))
    {
      image->Y.num_rows = image->U.num_rows = image->V.num_rows = num_rows;
      image->Y.num_cols = image->U.num_cols = image->V.num_cols = num_cols;
    }
  else
    {
      image->Y.num_rows = num_rows;
      image->U.num_rows = image->V.num_rows = 0;
      image->Y.num_cols = num_cols;
      image->U.num_cols = image->V.num_cols = 0;
    }
  
  return(0);
}


int QccIMGImageSetSizeYUV(QccIMGImage *image,
                          int num_rows_Y, int num_cols_Y,
                          int num_rows_U, int num_cols_U,
                          int num_rows_V, int num_cols_V)
{
  if (image == NULL)
    return(0);
  
  image->Y.num_rows = num_rows_Y;
  image->Y.num_cols = num_cols_Y;
  
  image->U.num_rows = num_rows_U;
  image->U.num_cols = num_cols_U;
  
  image->V.num_rows = num_rows_V;
  image->V.num_cols = num_cols_V;
  
  return(0);
}


int QccIMGImageAlloc(QccIMGImage *image)
{
  if (image == NULL)
    return(0);
  
  if (QccIMGImageComponentAlloc(&(image->Y)))
    {
      printf("(QccIMGImageAlloc): Error calling QccIMGImageComponentAlloc()");
      goto Error;
    }
  if (QccIMGImageComponentAlloc(&(image->U)))
    {
      printf("(QccIMGImageAlloc): Error calling QccIMGImageComponentAlloc()");
      goto Error;
    }
  if (QccIMGImageComponentAlloc(&(image->V)))
    {
      printf("(QccIMGImageAlloc): Error calling QccIMGImageComponentAlloc()");
      goto Error;
    }
  
  return(0);
  
 Error:
  QccIMGImageComponentFree(&(image->Y));
  QccIMGImageComponentFree(&(image->U));
  QccIMGImageComponentFree(&(image->V));
  return(1);
}


void QccIMGImageFree(QccIMGImage *image)
{
  if (image == NULL)
    return;
  
  QccIMGImageComponentFree(&(image->Y));
  QccIMGImageComponentFree(&(image->U));
  QccIMGImageComponentFree(&(image->V));
}


static int QccIMGImageSetMin(QccIMGImage *image)
{
  if (image == NULL)
    return(0);
  
  QccIMGImageComponentSetMin(&(image->Y));
  QccIMGImageComponentSetMin(&(image->U));
  QccIMGImageComponentSetMin(&(image->V));
  
  return(0);
}


static int QccIMGImageSetMax(QccIMGImage *image)
{
  if (image == NULL)
    return(0);
  
  QccIMGImageComponentSetMax(&(image->Y));
  QccIMGImageComponentSetMax(&(image->U));
  QccIMGImageComponentSetMax(&(image->V));
  
  return(0);
}


int QccIMGImageSetMaxMin(QccIMGImage *image)
{
  if (image == NULL)
    return(0);
  
  if (QccIMGImageSetMax(image))
    {
      printf("(QccIMGImageSetMaxMin): Error calling QccIMGImageSetMax()");
      return(1);
    }
  
  if (QccIMGImageSetMin(image))
    {
      printf("(QccIMGImageSetMaxMin): Error calling QccIMGImageSetMin()");
      return(1);
    }
  
  return(0);
}


int QccIMGImageColor(const QccIMGImage *image)
{
  return(image->image_type == QCCIMGTYPE_PPM);
}


int QccIMGImageDetermineType(QccIMGImage *image)
{
  QccString magic_num;
  QccString extension;
  
  if (image == NULL)
    return(0);
  
  if (image->image_type != QCCIMGTYPE_UNKNOWN)
    return(0);
  
  if (QccFileExists(image->filename))
    {
      if (QccFileGetMagicNumber(image->filename, magic_num))
        {
          printf("(QccIMGImageDetermineType): Error calling QccFileGetMagicNumber()");
          goto Return;
        }
      
      if (!strncmp(magic_num, "ICP", QCCSTRINGLEN))
        {
          image->image_type = QCCIMGTYPE_ICP;
          return(0);
        }
      
    }
  
  if (QccFileGetExtension(image->filename, extension))
    goto Error;
  
  if (!strncmp(extension, "icp", QCCSTRINGLEN))
    {
      image->image_type = QCCIMGTYPE_ICP;
      return(0);
    }
  if (!strncmp(extension, "pbm", QCCSTRINGLEN))
    {
      image->image_type = QCCIMGTYPE_PBM;
      return(0);
    }
  if (!strncmp(extension, "pgm", QCCSTRINGLEN))
    {
      image->image_type = QCCIMGTYPE_PGM;
      return(0);
    }
  if (!strncmp(extension, "ppm", QCCSTRINGLEN))
    {
      image->image_type = QCCIMGTYPE_PPM;
      return(0);
    }
  
 Error:
  image->image_type = QCCIMGTYPE_UNKNOWN;
  printf("(QccIMGImageDetermineType): Unknown image type");
 Return:
  return(1);
}


int QccIMGImageRead(QccIMGImage *image)
{
  if (image == NULL)
    return(0);
  
  image->image_type = QCCIMGTYPE_UNKNOWN;
  
  if (QccIMGImageDetermineType(image))
    {
      printf("(QccIMGImageRead): Error calling QccIMGImageDetermineType()");
      return(1);
    }
  
  switch (image->image_type)
    {
    case QCCIMGTYPE_ICP:
      QccStringCopy(image->Y.filename, image->filename);
      if (QccIMGImageComponentRead(&(image->Y)))
        {
          printf("(QccIMGImageRead): Error calling QccIMGImageComponentRead()");
          return(1);
        }
      break;
/*      
    case QCCIMGTYPE_PBM:
    case QCCIMGTYPE_PGM:
    case QCCIMGTYPE_PPM:
      if (QccIMGImagePNMRead(image))
        {
          printf("(QccIMGImageRead): Error calling QccIMGImagePNMRead()");
          return(1);
        }
      break;
   */   
    default:
      printf("(QccIMGImageRead): Unknown image type");
      return(1);
    }
  
  if (QccIMGImageSetMaxMin(image))
    {
      printf("(QccIMGImageRead): Error calling QccIMGImageSetMaxMin()");
      return(1);
    }
  
  return(0);
}



int QccIMGImageWrite(QccIMGImage *image)
{
  if (image == NULL)
    return(0);
  
  if (QccIMGImageDetermineType(image))
    {
      printf("(QccIMGImageWrite): Error calling QccIMGImageDetermineType()");
      return(1);
    }
  
  switch (image->image_type)
    {
    case QCCIMGTYPE_ICP:
      QccStringCopy(image->Y.filename, image->filename);
      if (QccIMGImageComponentWrite(&(image->Y)))
        {
          printf("(QccIMGImageWrite): Error calling QccIMGImageComponentWrite()");
          return(1);
        }
      break;
  /*    
    case QCCIMGTYPE_PBM:
    case QCCIMGTYPE_PGM:
    case QCCIMGTYPE_PPM:
      if (QccIMGImagePNMWrite(image))
        {
          printf("(QccIMGImageWrite): Error calling QccIMGImagePNMWrite()");
          return(1);
        }
      break;
     */ 
    default:
      printf("(QccIMGImageWrite): Unknown image type");
      return(1);
    }
  
  return(0);
}


/*  motion_estimation.c  */
#define QCCVID_ME_FULLPIXEL 0
#define QCCVID_ME_HALFPIXEL 1
#define QCCVID_ME_QUARTERPIXEL 2
#define QCCVID_ME_EIGHTHPIXEL 3



int QccVIDMotionEstimationExtractBlock(const QccIMGImageComponent *image,
                                       double row,
                                       double col,
                                       QccMatrix block,
                                       int block_size,
                                       int subpixel_accuracy)
{
  int block_row, block_col;
  int step;
  int row2, col2;
  int extraction_row, extraction_col;

  if (image == NULL)
    return(0);
  if (image->image == NULL)
    return(0);
  if (block == NULL)
    return(0);

  switch (subpixel_accuracy)
    {
    case QCCVID_ME_FULLPIXEL:
      step = 1;
      row2 = (int)row;
      col2 = (int)col;
      break;
    case QCCVID_ME_HALFPIXEL:
      step = 2;
      row2 = (int)(row * 2);
      col2 = (int)(col * 2);
      break;
    case QCCVID_ME_QUARTERPIXEL:
      step = 4;
      row2 = (int)(row * 4);
      col2 = (int)(col * 4);
      break;
    case QCCVID_ME_EIGHTHPIXEL:
      step = 8;
      row2 = (int)(row * 8);
      col2 = (int)(col * 8);
      break;
    default:
      printf("(QccVIDMotionEstimationExtractBlock): Unrecognized subpixel accuracy (%d)",
                         subpixel_accuracy);
      return(1);
    }

  for (block_row = 0; block_row < block_size; block_row++)
    for (block_col = 0; block_col < block_size; block_col++)
      {
        extraction_row = block_row * step + row2;
        extraction_col = block_col * step + col2;

        if (extraction_row < 0)
          extraction_row = 0;
        else
          if (extraction_row >= image->num_rows)
            extraction_row = image->num_rows - 1;

        if (extraction_col < 0)
          extraction_col = 0;
        else
          if (extraction_col >= image->num_cols)
            extraction_col = image->num_cols - 1;

        block[block_row][block_col] =
          image->image[extraction_row][extraction_col];
      }

  return(0);
}


int QccVIDMotionEstimationInsertBlock(QccIMGImageComponent *image,
                                      double row,
                                      double col,
                                      const QccMatrix block,
                                      int block_size,
                                      int subpixel_accuracy)
{
  int block_row, block_col;
  int step;
  int row2, col2;
  int insertion_row, insertion_col;

  if (image == NULL)
    return(0);
  if (image->image == NULL)
    return(0);
  if (block == NULL)
    return(0);

  switch (subpixel_accuracy)
    {
    case QCCVID_ME_FULLPIXEL:
      step = 1;
      row2 = (int)row;
      col2 = (int)col;
      break;
    case QCCVID_ME_HALFPIXEL:
      step = 2;
      row2 = (int)(row * 2);
      col2 = (int)(col * 2);
      break;
    case QCCVID_ME_QUARTERPIXEL:
      step = 4;
      row2 = (int)(row * 4);
      col2 = (int)(col * 4);
      break;
    case QCCVID_ME_EIGHTHPIXEL:
      step = 8;
      row2 = (int)(row * 8);
      col2 = (int)(col * 8);
      break;
    default:
      printf("(QccVIDMotionEstimationInsertBlock): Unrecognized subpixel accuracy (%d)",
                         subpixel_accuracy);
      return(1);
    }

  for (block_row = 0; block_row < block_size; block_row++)
    for (block_col = 0; block_col < block_size; block_col++)
      {
        insertion_row = block_row * step + row2;
        insertion_col = block_col * step + col2;

        if ((insertion_row >= 0) &&
            (insertion_row < image->num_rows) &&
            (insertion_col >= 0) &&
            (insertion_col < image->num_cols))
          image->image[insertion_row][insertion_col] =
            block[block_row][block_col];
      }

  return(0);
}


double QccVIDMotionEstimationMAE(QccMatrix current_block,
                                 QccMatrix reference_block,
                                 QccMatrix weights,
                                 int block_size)
{
  double mae = 0.0;
  int block_row, block_col;

  if (current_block == NULL)
    return(0.0);
  if (reference_block == NULL)
    return(0.0);

  if (weights == NULL)
    for (block_row = 0; block_row < block_size; block_row++)
      for (block_col = 0; block_col < block_size; block_col++)
        mae +=
          fabs(current_block[block_row][block_col] -
               reference_block[block_row][block_col]);
  else
    for (block_row = 0; block_row < block_size; block_row++)
      for (block_col = 0; block_col < block_size; block_col++)
        mae +=
          weights[block_row][block_col] *
          fabs(current_block[block_row][block_col] -
               reference_block[block_row][block_col]);

  return(mae / block_size / block_size);
}


static int QccVIDMotionEstimationWindowSearch(QccMatrix current_block,
                                              QccMatrix reference_block,
                                              int block_size,
                                              const QccIMGImageComponent
                                              *reference_frame,
                                              double search_row,
                                              double search_col,
                                              double window_size,
                                              double search_step,
                                              int subpixel_accuracy,
                                              double *motion_vector_horizontal,
                                              double *motion_vector_vertical)
{
  double u, v;
  double reference_frame_row, reference_frame_col;
  double current_mae;
  double min_mae = MAXDOUBLE;

  for (v = -window_size; v <= window_size; v += search_step)
    for (u = -window_size; u <= window_size; u += search_step)
      {
        reference_frame_row = search_row + v;
        reference_frame_col = search_col + u;
        
        if (QccVIDMotionEstimationExtractBlock(reference_frame,
                                               reference_frame_row,
                                               reference_frame_col,
                                               reference_block,
                                               block_size,
                                               subpixel_accuracy))
          {
            printf("(QccVIDMotionEstimationWindowSearch): Error calling QccVIDMotionEstimationExtractBlock()");
            return(1);
          }
        
        current_mae = QccVIDMotionEstimationMAE(current_block,
                                                reference_block,
                                                NULL,
                                                block_size);
        
        if (current_mae < min_mae)
          {
            min_mae = current_mae;
            *motion_vector_horizontal = u;
            *motion_vector_vertical = v;
          }
      }

  return(0);
}


int QccVIDMotionEstimationFullSearch(const QccIMGImageComponent *current_frame,
                                     const QccIMGImageComponent
                                     *reference_frame,
                                     QccIMGImageComponent
                                     *motion_vectors_horizontal,
                                     QccIMGImageComponent
                                     *motion_vectors_vertical,
                                     int block_size,
                                     int window_size,
                                     int subpixel_accuracy)
{
  int return_value;
  QccMatrix current_block = NULL;
  QccMatrix reference_block = NULL;
  int block_row, block_col;
  double u = 0;
  double v = 0;
  double current_search_step;
  double final_search_step;
  int mv_row, mv_col;
  double search_row, search_col;
  double current_window_size;

  if (current_frame == NULL)
    return(0);
  if (reference_frame == NULL)
    return(0);
  if (motion_vectors_horizontal == NULL)
    return(0);
  if (motion_vectors_vertical == NULL)
    return(0);

  if (current_frame->image == NULL)
    return(0);
  if (reference_frame->image == NULL)
    return(0);
  if (motion_vectors_horizontal->image == NULL)
    return(0);
  if (motion_vectors_vertical->image == NULL)
    return(0);

  if ((motion_vectors_horizontal->num_rows !=
       current_frame->num_rows / block_size) ||
      (motion_vectors_horizontal->num_cols !=
       current_frame->num_cols / block_size))
    {
      printf("(QccVIDMotionEstimationFullSearch): Motion-vector field is inconsistent with current-frame size");
      goto Error;
    }

  if ((motion_vectors_vertical->num_rows !=
       current_frame->num_rows / block_size) ||
      (motion_vectors_vertical->num_cols !=
       current_frame->num_cols / block_size))
    {
      printf("(QccVIDMotionEstimationFullSearch): Motion-vector field is inconsistent with current-frame size");
      goto Error;
    }

  switch (subpixel_accuracy)
    {
    case QCCVID_ME_FULLPIXEL:
      if ((reference_frame->num_rows != current_frame->num_rows) ||
          (reference_frame->num_cols != current_frame->num_cols))
        {
          printf("(QccVIDMotionEstimationFullSearch): Reference-frame size is inconsistent with current-frame size for full-pixel motion estimation");
          goto Error;
        }
      final_search_step = 1.0;
      break;
    case QCCVID_ME_HALFPIXEL:
      if ((reference_frame->num_rows != 2 * current_frame->num_rows) ||
          (reference_frame->num_cols != 2 * current_frame->num_cols))
        {
          printf("(QccVIDMotionEstimationFullSearch): Reference-frame size is inconsistent with current-frame size for half-pixel motion estimation");
          goto Error;
        }
      final_search_step = 0.5;
      break;
    case QCCVID_ME_QUARTERPIXEL:
      if ((reference_frame->num_rows != 4 * current_frame->num_rows) ||
          (reference_frame->num_cols != 4 * current_frame->num_cols))
        {
          printf("(QccVIDMotionEstimationFullSearch): Reference-frame size is inconsistent with current-frame size for quarter-pixel motion estimation");
          goto Error;
        }
      final_search_step = 0.25;
      break;
    case QCCVID_ME_EIGHTHPIXEL:
      if ((reference_frame->num_rows != 8 * current_frame->num_rows) ||
          (reference_frame->num_cols != 8 * current_frame->num_cols))
        {
          printf("(QccVIDMotionEstimationFullSearch): Reference-frame size is inconsistent with current-frame size for eighth-pixel motion estimation");
          goto Error;
        }
      final_search_step = 0.125;
      break;
    default:
      printf("(QccVIDMotionEstimationFullSearch): Unrecognized motion-estimation accuracy (%d)",
                         subpixel_accuracy);
      goto Error;
    }

  if ((current_block =
       QccMatrixAlloc(block_size, block_size)) == NULL)
    {
      printf("(QccVIDMotionEstimationFullSearch): Error calling QccMatrixAlloc()");
      goto Error;
    }
  if ((reference_block =
       QccMatrixAlloc(block_size, block_size)) == NULL)
    {
      printf("(QccVIDMotionEstimationFullSearch): Error calling QccMatrixAlloc()");
      goto Error;
    }

  for (block_row = 0;
       block_row < current_frame->num_rows; block_row += block_size)
    for (block_col = 0;
         block_col < current_frame->num_cols; block_col += block_size)
      {
        if (QccVIDMotionEstimationExtractBlock(current_frame,
                                               (double)block_row,
                                               (double)block_col,
                                               current_block,
                                               block_size,
                                               QCCVID_ME_FULLPIXEL))
          {
            printf("(QccVIDMotionEstimationFullSearch): Error calling QccVIDMotionEstimationExtractBlock()");
            goto Error;
          }

        mv_row = block_row / block_size;
        mv_col = block_col / block_size;

        motion_vectors_horizontal->image[mv_row][mv_col] = 0.0;
        motion_vectors_vertical->image[mv_row][mv_col] = 0.0;

        for (current_search_step = 1.0;
             current_search_step >= final_search_step;
             current_search_step /= 2)
          {
            search_row = block_row + 
              motion_vectors_vertical->image[mv_row][mv_col];
            search_col = block_col + 
              motion_vectors_horizontal->image[mv_row][mv_col];

            current_window_size =
              (current_search_step == 1.0) ?
              (double)window_size : current_search_step;

            if (QccVIDMotionEstimationWindowSearch(current_block,
                                                   reference_block,
                                                   block_size,
                                                   reference_frame,
                                                   search_row,
                                                   search_col,
                                                   current_window_size,
                                                   current_search_step,
                                                   subpixel_accuracy,
                                                   &u,
                                                   &v))
              {
                printf("(QccVIDMotionEstimationFullSearch): Error calling QccVIDMotionEstimationWindowSearch()");
                goto Error;
              }

            motion_vectors_horizontal->image[mv_row][mv_col] += u;
            motion_vectors_vertical->image[mv_row][mv_col] += v;
          }
      }
      
  return_value = 0;
  goto Return;
 Error:
  return_value = 1;
 Return:
  QccMatrixFree(current_block, block_size);
  QccMatrixFree(reference_block, block_size);
  return(return_value);
}


int QccVIDMotionEstimationCalcReferenceFrameSize(int num_rows,
                                                 int num_cols,
                                                 int *reference_num_rows,
                                                 int *reference_num_cols,
                                                 int subpixel_accuracy)
{
  switch (subpixel_accuracy)
    {
    case QCCVID_ME_FULLPIXEL:
      *reference_num_rows = num_rows;
      *reference_num_cols = num_cols;
      break;
    case QCCVID_ME_HALFPIXEL:
      *reference_num_rows = num_rows * 2;
      *reference_num_cols = num_cols * 2;
      break;
    case QCCVID_ME_QUARTERPIXEL:
      *reference_num_rows = num_rows * 4;
      *reference_num_cols = num_cols * 4;
      break;
    case QCCVID_ME_EIGHTHPIXEL:
      *reference_num_rows = num_rows * 8;
      *reference_num_cols = num_cols * 8;
      break;
    default:
      printf("(QccVIDMotionEstimationCalcReferenceFrameSize): Unrecognized subpixel accuracy (%d)",
                         subpixel_accuracy);
      return(1);
    }

  return(0);
}


int QccVIDMotionEstimationCreateReferenceFrame(const QccIMGImageComponent
                                               *current_frame,
                                               QccIMGImageComponent
                                               *reference_frame,
                                               int subpixel_accuracy,
                                               const QccFilter *filter1,
                                               const QccFilter *filter2,
                                               const QccFilter *filter3)
{
  int return_value;
  QccIMGImageComponent reference_frame2;
  QccIMGImageComponent reference_frame3;

  QccIMGImageComponentInitialize(&reference_frame2);
  QccIMGImageComponentInitialize(&reference_frame3);

  if (current_frame == NULL)
    return(0);
  if (reference_frame == NULL)
    return(0);

  switch (subpixel_accuracy)
    {
    case QCCVID_ME_FULLPIXEL:
      if ((reference_frame->num_rows != current_frame->num_rows) ||
          (reference_frame->num_cols != current_frame->num_cols))
        {
          printf("(QccVIDMotionEstimationCreateReferenceFrame): Reference-frame size is inconsistent with current-frame size for full-pixel motion estimation");
          goto Error;
        }
      if (QccIMGImageComponentCopy(reference_frame,
                                   current_frame))
        {
          printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentCopy()");
          goto Error;
        }
      break;
    case QCCVID_ME_HALFPIXEL:
      if ((reference_frame->num_rows != 2 * current_frame->num_rows) ||
          (reference_frame->num_cols != 2 * current_frame->num_cols))
        {
          printf("(QccVIDMotionEstimationCreateReferenceFrame): Reference-frame size is inconsistent with current-frame size for half-pixel motion estimation");
          goto Error;
        }
      if (filter1 == NULL)
        {
          if (QccIMGImageComponentInterpolateBilinear(current_frame,
                                                      reference_frame))
            {
              printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentInterpolateBilinear()");
              goto Error;
            }
        }
      else
        if (QccIMGImageComponentInterpolateFilter(current_frame,
                                                  reference_frame,
                                                  filter1))
          {
            printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentInterpolateFilter()");
            goto Error;
          }
      break;
    case QCCVID_ME_QUARTERPIXEL:
      if ((reference_frame->num_rows != 4 * current_frame->num_rows) ||
          (reference_frame->num_cols != 4 * current_frame->num_cols))
        {
          printf("(QccVIDMotionEstimationCreateReferenceFrame): Reference-frame size is inconsistent with current-frame size for quarter-pixel motion estimation");
          goto Error;
        }
      reference_frame2.num_rows = 2 * current_frame->num_rows;
      reference_frame2.num_cols = 2 * current_frame->num_cols;
      if (QccIMGImageComponentAlloc(&reference_frame2))
        {
          printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentAlloc()");
          goto Error;
        }
      if (filter1 == NULL)
        {
          if (QccIMGImageComponentInterpolateBilinear(current_frame,
                                                      &reference_frame2))
            {
              printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentInterpolateBilinear()");
              goto Error;
            }
        }
      else
        if (QccIMGImageComponentInterpolateFilter(current_frame,
                                                  &reference_frame2,
                                                  filter1))
          {
            printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentInterpolateFilter()");
            goto Error;
          }
      if (filter2 == NULL)
        {
          if (QccIMGImageComponentInterpolateBilinear(&reference_frame2,
                                                      reference_frame))
            {
              printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentInterpolateBilinear()");
              goto Error;
            }
        }
      else
        if (QccIMGImageComponentInterpolateFilter(&reference_frame2,
                                                  reference_frame,
                                                  filter2))
          {
            printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentInterpolateFilter()");
            goto Error;
          }
      break;
    case QCCVID_ME_EIGHTHPIXEL:
      if ((reference_frame->num_rows != 8 * current_frame->num_rows) ||
          (reference_frame->num_cols != 8 * current_frame->num_cols))
        {
          printf("(QccVIDMotionEstimationCreateReferenceFrame): Reference-frame size is inconsistent with current-frame size for eighth-pixel motion estimation");
          goto Error;
        }
      reference_frame2.num_rows = 2 * current_frame->num_rows;
      reference_frame2.num_cols = 2 * current_frame->num_cols;
      if (QccIMGImageComponentAlloc(&reference_frame2))
        {
          printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentAlloc()");
          goto Error;
        }
      reference_frame3.num_rows = 4 * current_frame->num_rows;
      reference_frame3.num_cols = 4 * current_frame->num_cols;
      if (QccIMGImageComponentAlloc(&reference_frame3))
        {
          printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentAlloc()");
          goto Error;
        }
      if (filter1 == NULL)
        {
          if (QccIMGImageComponentInterpolateBilinear(current_frame,
                                                      &reference_frame2))
            {
              printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentInterpolateBilinear()");
              goto Error;
            }
        }
      else
        if (QccIMGImageComponentInterpolateFilter(current_frame,
                                                  &reference_frame2,
                                                  filter1))
          {
            printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentInterpolateFilter()");
            goto Error;
          }
      if (filter2 == NULL)
        {
          if (QccIMGImageComponentInterpolateBilinear(&reference_frame2,
                                                      &reference_frame3))
            {
              printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentInterpolateBilinear()");
              goto Error;
            }
        }
      else
        if (QccIMGImageComponentInterpolateFilter(&reference_frame2,
                                                  &reference_frame3,
                                                  filter2))
          {
            printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentInterpolateFilter()");
            goto Error;
          }
      if (filter3 == NULL)
        {
          if (QccIMGImageComponentInterpolateBilinear(&reference_frame3,
                                                      reference_frame))
            {
              printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentInterpolateBilinear()");
              goto Error;
            }
        }
      else
        if (QccIMGImageComponentInterpolateFilter(&reference_frame3,
                                                  reference_frame,
                                                  filter3))
          {
            printf("(QccVIDMotionEstimationCreateReferenceFrame): Error calling QccIMGImageComponentInterpolateFilter()");
            goto Error;
          }
      break;
    default:
      printf("(QccVIDMotionEstimationCreateReferenceFrame): Unrecognized motion-estimation accuracy (%d)",
                         subpixel_accuracy);
      goto Error;
    }

  return_value = 0;
  goto Return;
 Error:
  return_value = 1;
 Return:
  QccIMGImageComponentFree(&reference_frame2);
  QccIMGImageComponentFree(&reference_frame3);
  return(return_value);
}


int QccVIDMotionEstimationCreateCompensatedFrame(QccIMGImageComponent
                                                 *motion_compensated_frame,
                                                 const QccIMGImageComponent
                                                 *reference_frame,
                                                 const QccIMGImageComponent
                                                 *motion_vectors_horizontal,
                                                 const QccIMGImageComponent
                                                 *motion_vectors_vertical,
                                                 int block_size,
                                                 int subpixel_accuracy)
{
  int return_value;
  QccMatrix reference_block = NULL;
  int block_row, block_col;
  int mv_row, mv_col;

  if (motion_compensated_frame == NULL)
    return(0);
  if (reference_frame == NULL)
    return(0);
  if (motion_vectors_horizontal == NULL)
    return(0);
  if (motion_vectors_vertical == NULL)
    return(0);

  if (motion_compensated_frame->image == NULL)
    return(0);
  if (reference_frame->image == NULL)
    return(0);
  if (motion_vectors_horizontal->image == NULL)
    return(0);
  if (motion_vectors_vertical->image == NULL)
    return(0);

  if ((motion_vectors_horizontal->num_rows !=
       motion_compensated_frame->num_rows / block_size) ||
      (motion_vectors_horizontal->num_cols !=
       motion_compensated_frame->num_cols / block_size))
    {
      printf("(QccVIDMotionEstimationCreateCompensatedFrame): Motion-vector field is inconsistent with current-frame size");
      goto Error;
    }

  if ((motion_vectors_vertical->num_rows !=
       motion_compensated_frame->num_rows / block_size) ||
      (motion_vectors_vertical->num_cols !=
       motion_compensated_frame->num_cols / block_size))
    {
      printf("(QccVIDMotionEstimationCreateCompensatedFrame): Motion-vector field is inconsistent with current-frame size");
      goto Error;
    }

  switch (subpixel_accuracy)
    {
    case QCCVID_ME_FULLPIXEL:
      if ((reference_frame->num_rows != motion_compensated_frame->num_rows) ||
          (reference_frame->num_cols != motion_compensated_frame->num_cols))
        {
          printf("(QccVIDMotionEstimationCreateCompensatedFrame): Reference-frame size is inconsistent with current-frame size for full-pixel motion estimation");
          goto Error;
        }
      break;
    case QCCVID_ME_HALFPIXEL:
      if ((reference_frame->num_rows !=
           2 * motion_compensated_frame->num_rows) ||
          (reference_frame->num_cols !=
           2 * motion_compensated_frame->num_cols))
        {
          printf("(QccVIDMotionEstimationCreateCompensatedFrame): Reference-frame size is inconsistent with current-frame size for half-pixel motion estimation");
          goto Error;
        }
      break;
    case QCCVID_ME_QUARTERPIXEL:
      if ((reference_frame->num_rows !=
           4 * motion_compensated_frame->num_rows) ||
          (reference_frame->num_cols !=
           4 * motion_compensated_frame->num_cols))
        {
          printf("(QccVIDMotionEstimationCreateCompensatedFrame): Reference-frame size is inconsistent with current-frame size for quarter-pixel motion estimation");
          goto Error;
        }
      break;
    case QCCVID_ME_EIGHTHPIXEL:
      if ((reference_frame->num_rows !=
           8 * motion_compensated_frame->num_rows) ||
          (reference_frame->num_cols !=
           8 * motion_compensated_frame->num_cols))
        {
          printf("(QccVIDMotionEstimationCreateCompensatedFrame): Reference-frame size is inconsistent with current-frame size for eighth-pixel motion estimation");
          goto Error;
        }
      break;
    default:
      printf("(QccVIDMotionEstimationCreateCompensatedFrame): Unrecognized motion-estimation accuracy (%d)",
                         subpixel_accuracy);
      goto Error;
    }

  if ((reference_block = QccMatrixAlloc(block_size, block_size)) == NULL)
    {
      printf("(QccVIDMotionEstimationCreateCompensatedFrame): Error calling QccMatrixAlloc()");
      goto Error;
    }

  for (block_row = 0;
       block_row < motion_compensated_frame->num_rows;
       block_row += block_size)
    for (block_col = 0;
         block_col < motion_compensated_frame->num_cols;
         block_col += block_size)
      {
        mv_row = block_row / block_size;
        mv_col = block_col / block_size;

        if (QccVIDMotionEstimationExtractBlock(reference_frame,
                                               (double)block_row +
                                               motion_vectors_vertical->image
                                               [mv_row][mv_col],
                                               (double)block_col +
                                               motion_vectors_horizontal->image
                                               [mv_row][mv_col],
                                               reference_block,
                                               block_size,
                                               subpixel_accuracy))
          {
            printf("(QccVIDMotionEstimationCreateCompensatedFrame): Error calling QccVIDMotionEstimationExtractBlock()");
            goto Error;
          }

        if (QccVIDMotionEstimationInsertBlock(motion_compensated_frame,
                                              (double)block_row,
                                              (double)block_col,
                                              reference_block,
                                              block_size,
                                              QCCVID_ME_FULLPIXEL))
          {
            printf("(QccVIDMotionEstimationCreateCompensatedFrame): Error calling QccVIDMotionEstimationInsertBlock()");
            goto Error;
          }
      }

  return_value = 0;
  goto Return;
 Error:
  return_value = 1;
 Return:
  QccMatrixFree(reference_block, block_size);
  return(return_value);
}


#endif /* MEMC_H */
