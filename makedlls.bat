REM
REM ## Makefile for generating the DLL files censre4d.dll and surv_leftcens_cov_v1.dll ##
REM


REM ########################################################################################

REM
REM Make the 32-bit windows DLLs - i386 architecture
REM

gcc -m32 -shared -c csubs.c -I "C:\Users\xxie\Documents\R\R-3.6.2\include"
gfortran -m32 -shared -fno-optimize-sibling-calls -c censre4d.f

gcc -m32 -shared -o censre4d_i386.dll censre4d.o csubs.o -lR -lRlapack -lRblas -L "C:\Users\xxie\Documents\R\R-3.6.2\bin\i386"


gfortran -m32 -shared -fno-optimize-sibling-calls -c censre_subs.f
gfortran -m32 -shared -fno-optimize-sibling-calls -c surv_left_cens_cov_estc.f
gfortran -m32 -shared -fno-optimize-sibling-calls -c surv_left_cens_cov_infa.f

gfortran -m32 -shared -o surv_leftcens_cov_v1_i386.dll surv_left_cens_cov_estc.o surv_left_cens_cov_infa.o censre_subs.o csubs.o -lR -lRlapack -lRblas -L "C:\Users\xxie\Documents\R\R-3.6.2\bin\i386"


del *.o


REM ########################################################################################

REM
REM Make the 64-bit windows DLLs - i386 architecture
REM

gcc -m64 -shared -c csubs.c -I "C:\Program Files\R\R-4.4.1\include"
gfortran -m64 -shared -fno-optimize-sibling-calls -c censre4d.f

gcc -m64 -shared -o censre4d_x64.dll censre4d.o csubs.o -lR -lRlapack -lRblas -L "C:\Program Files\R\R-4.4.1\bin\x64"


gfortran -m64 -shared -fno-optimize-sibling-calls -c censre_subs.f
gfortran -m64 -shared -fno-optimize-sibling-calls -c surv_left_cens_cov_estc.f
gfortran -m64 -shared -fno-optimize-sibling-calls -c surv_left_cens_cov_infa.f

gfortran -m64 -shared -o surv_leftcens_cov_v1_x64.dll surv_left_cens_cov_estc.o surv_left_cens_cov_infa.o censre_subs.o csubs.o -lR -lRlapack -lRblas -L "C:\Program Files\R\R-4.4.1\bin\x64"


del *.o

REM ########################################################################################
