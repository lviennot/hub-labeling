
build: _build
	cd $<; make -j

_build:
	mkdir -p $@
	cd $@; ln -sf ../test_data; cmake ..

unit: build
	_build/unit

test: _data/iceland-latest.osm-d.gr.gz build
	@echo "\nThis test takes less than a minute on my laptop...\n"
	gunzip -c $< | awk '{print($$2,$$3,$$4);}' | _build/hl_trans test -


REPO:=https://files.inria.fr/gang/graphs/osm2015road/europe/

_data/%:
	mkdir -p _data
	curl -o $@ $(REPO)/$*


clean:
	rm -f *.o src/*~ *~
	rm -fr *.dSYM _*

%.force:
	rm -f $*; make $*

