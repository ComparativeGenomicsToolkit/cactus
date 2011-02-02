all :
	cd src && make all

clean : clean.src
	rm -rf lib/*.h bin/*.dSYM 
	
clean.src :
	cd src && make clean

