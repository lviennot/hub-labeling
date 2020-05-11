HDRS:=$(wildcard src/*.hh)
CTRS:=$(wildcard src/*.cc)
UNTS:=$(wildcard src/*_unit.cc)
LIBS:=

main: hl_trans.o
	@echo "Try it with 'make test'"

test: _data/iceland-latest.osm-d.gr.gz hl_trans.o
	@echo "\nThis test takes less than a minute on my laptop...\n"
	gunzip -c $< | awk '{print($$2,$$3,$$4);}' | ./hl_trans.o test -

unit: src/unit.cc $(HDRS) $(UNTS)
	g++ -std=c++11 -g -rdynamic -pthread $(LIBS) -o $@ $(UNTS) $< 
	./unit

%.o: src/%.cc $(HDRS)
	g++ -std=c++11 -O3 -pthread $(LIBS) -o $@ $<


REPO:=https://files.inria.fr/gang/graphs/osm2015road/europe/

_data/%:
	mkdir -p _data
	curl -o $@ $(REPO)/$*

clean:
	rm -fr *~ */*~ _* *.o
