CPPPATH= -I '/home/gandhali/cantera-May26/canteraBuild/include/' -I '/home/gandhali/anaconda3/envs/ct-build/include' -I '/home/gandhali/anaconda3/envs/ct-build/include/'
LIBS= -L'/home/gandhali/cantera-May26/canteraBuild/lib' -L'/home/gandhali/anaconda3/envs/ct-build/lib' -lboost_filesystem -lcantera -lsundials_cvodes -lsundials_ida -lsundials_nvecserial
LINKFLAGS=['/nologo']

all: full 

full: 
	g++ $(LIKFLAGS) $(CPPPATH) -o boundaryLayer main.cpp solver.cpp BoundaryLayer.cpp -pthread -O3  $(LIBS) 

clean: 
	rm output 
