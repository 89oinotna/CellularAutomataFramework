FF_ROOT	= -I/home/kkk/fastflow
LDFLAGS	=  -pthread -lX11 -lpng -O3 -finline-functions
CXX = g++-10 
IMG = -DWIMG
TARGETS = mine sequential ff_parfor ff_farm

$(TARGETS): %: %.cpp
	$(CXX) $<  $(LDFLAGS) $(FF_ROOT) -o $@

all: $(TARGETS) sequentialw minew ff_parforw ff_farmw

sequentialw: sequential.cpp
	$(CXX) sequential.cpp $(LDFLAGS) $(IMG) -o sequential_write

minew: mine.cpp
	$(CXX) mine.cpp  $(LDFLAGS) $(FF_ROOT) $(IMG) -o mine_write

ff_parforw: ff_parfor.cpp
	$(CXX) ff_parfor.cpp  $(LDFLAGS) $(FF_ROOT) $(IMG) -o ff_parfor_write

ff_farmw:ff_farm.cpp
	$(CXX) ff_farm.cpp  $(LDFLAGS) $(FF_ROOT) $(IMG) -o ff_farm_write

	
