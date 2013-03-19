
all:	buddhabrot

buddhabrot:	Buddhabrot.cpp
	g++ -g3 -ggdb -O0 -DDEBUG Buddhabrot.cpp -o buddhabrot

