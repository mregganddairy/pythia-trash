tutorial%: tutorial%.cc
	g++ -I/usr/share/pythia8/include `root-config --cflags` $@.cc -o $@ -lpythia8 -L/usr/share/pythia8/lib `root-config --glibs`
