SRC_FILES := $(shell find -name *.cpp)
LDFLAGS := -lX11 -lGL -lpthread -lpng -lstdc++fs
CXXFLAGS := -Wall -Wextra -pedantic
DEBUGFLAGS := --debug -DDEBUG
RELEASEFLAGS := -Ofast -march=native -mtune=native

debug: $(SRC_FILES)
	g++ $(CXXFLAGS) $(DEBUGFLAGS) -o $@ $^ $(LDFLAGS)

release: $(SRC_FILES)
	g++ $(CXXFLAGS) $(RELEASEFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -rf release debug
