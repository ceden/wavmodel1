

include ../make.inc
code   = ../for_src


all: 
	echo "do nothing"

%.o: %.f90
	cp $< ${code}/config.f90
	(cd ${code};make all )

%: %.o
	mv ${code}/model.x ./${@}.x

clean: 
	(cd ${code}; make clean)


