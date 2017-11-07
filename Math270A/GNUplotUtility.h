//######################################################################
//   Class GNUplotUtility
//######################################################################
/*
 Class GNUplotUtility : A class whose member functions output 2D
 data into a gnuplot readable files.
 
 To view a wireframe surface representation of the data,
 start gnuplot and use a command such as
 
 splot "datafileName" w l
 
 To view a filled surface representation of the data, start
 gnuplot plot and use the the commands
 
 set pm3d
 splot "datafileName" w l
 
 (the set pm3d command need not be repeated for subsequent plots).
 
 */
// Math 270A
// Oct. 24, 2017
//######################################################################
//
//

// MS compilers generate warnings if fopen is used instead of fopen_s (a Microsoft specific language
// extension, so this macro implements the appropriate equivalent to fopen that MS wants when
// MS compilers are being used. In both versions, the
// macro returns a non-zero value if the open fails (e.g. a non-zero error code).
//

#ifndef _MSC_VER
#define OPENFILE(dataFile,filename,mode) ((dataFile = fopen(filename,  mode)) == NULL)
#else
#define OPENFILE(dataFile,fileName,mode) ((fopen_s(&dataFile,fileName, mode)) != 0)
#endif


#include "GridFun2D.h"

#ifndef _GNUplotUtility_
#define _GNUplotUtility_

class GNUplotUtility
{
public:
    
    static void output(GridFun2D& F, const string& fileName)
    {
        //
        //  Open and then write data the the gnuplot output file
        //
        FILE* dataFile;
        
        if(OPENFILE(dataFile,fileName.c_str(),"w"))  // Using OPENFILE macro defined at top of file.
        {
            printf("The file %s could not be opened\n",fileName.c_str());
            return;
        }
        
        
        double x; double y;
        for(long i = 0; i <= F.xPanel; i++)
        {
            x = F.xMin + i*F.hx;
            for(long j = 0; j <= F.yPanel; j++)
            {
                y =  F.yMin + j*F.hy;
                fprintf(dataFile,"%-10.5e %-10.5e %-10.5e \n",x,y,F.values(i,j));
            }
            fprintf(dataFile,"\n");
        }
        
        fclose(dataFile);
    }
    
};
#endif /* _GNUplotUtility_ */

