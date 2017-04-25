build: makedirs roxygen
	R CMD build KEADrugResponse
	mv KEADrugResponse_*.tar.gz build/

makedirs:
	mkdir -p build

clean:
	rm -rf build

roxygen:
	R CMD BATCH --no-save  pkgsrc/roxygenize.R /dev/stdout

vignette: makedirs
	cp -rf KEADrugResponse/vignettes build
	cp pkgsrc/vignette/Makefile build/vignettes
	make -C build/vignettes build

manual: makedirs
	R CMD Rd2pdf --no-preview KEADrugResponse
	mv -f KEADrugResponse.pdf build/

