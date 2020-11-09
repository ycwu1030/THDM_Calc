
LT = /Users/ycwu/Library/Mathematica/Applications/LoopTools/x86_64-Darwin/lib
#LT = /Users/mac/work/LoopTools/x86_64-Darwin/lib

SRCDIR := src
INCDIR := include
OBJDIR := obj
CXX = $(LT)/../bin/f++
FFLAG = -I$(LT)/../include -I$(INCDIR) -std=c++11
FLIBS = -L$(LT) -looptools

muSTULIBPATH := ../HiggsSignalStrength_STU
muSTUFLAGS := -I$(muSTULIBPATH)/
muSTULIBS := -L$(muSTULIBPATH)/ -lHiggsSignalStrengthSTU -Wl,-rpath,$(muSTULIBPATH)

FFLAG += $(muSTUFLAGS)
FLIBS += $(muSTULIBS)

SRC = $(wildcard $(SRCDIR)/*.cpp)
OBJ = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRC))

all: THDM_SCAN.x #THDMHCoupCalc.x THDMHCoupCalc_NC.x THDMHCoupCalc_Check.x

%.x:%.cpp $(OBJ)
	$(CXX) $(FFLAG) -o $@ $< $(OBJ) $(FLIBS)

THDMHCoupCalc.x: THDMHCoupCalc.cpp $(OBJ)
	$(CXX) $(FFLAG) -o $@ $< $(OBJ) $(FLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(FFLAG) -c $< -o $@

.PHONY: clean

clean:
	rm -f *.x
	rm -f $(OBJDIR)/*.o

