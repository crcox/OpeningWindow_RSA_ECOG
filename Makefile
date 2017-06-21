export MCC=$(MATLABDIR)/bin/mcc
export MEX=$(MATLABDIR)/bin/mex
export MFLAGS=-m -R -singleCompThread -R -nodisplay -R -nojvm
TOP := $(shell pwd)
SRCTAR=source_code.tar.gz
SRC=src
DEP=dependencies
DATA=data
GLMNET=$(DEP)/glmnet
INCL= -N -p ${MATLABDIR}/toolbox/stats/stats -I $(SRC) -I $(GLMNET)
.PHONEY: all clean-all clean-postbuild glmnet sdist

all: setup glmnet OpeningWindow_RSA_ECOG clean-postbuild

setup: $(SRC) $(DEP) $(DATA)
$(SRC) $(DEP) $(DATA):
	tar xzvf $(SRCTAR)

glmnet: $(GLMNET)/glmnetMex.mexa64
$(GLMNET)/glmnetMex.mexa64: $(GLMNET)/glmnetMex.F $(GLMNET)/GLMnet.f
	$(MEX) -fortran -outdir $(GLMNET) $^

OpeningWindow_RSA_ECOG: $(SRC)/OpeningWindow_RSA_ECOG.m
	$(MCC) -v $(MFLAGS) $(INCL) -a $(DATA) -o $@ $<

clean-postbuild:
	rm *.dmr
	rm mccExcludedFiles.log
	rm readme.txt

sdist:
	tar czhf $(SRCTAR) src dependencies data

clean-all:
	-rm OpeningWindow_RSA_ECOG
	-rm $(GLMNET)/glmnetMEX.mexa64
