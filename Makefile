SRC_FILES := $(shell find -name *.cpp)
LDFLAGS := -lX11 -lGL -lpthread -lpng -lstdc++fs
CXXFLAGS := -Wall -Wextra -pedantic
DEBUGFLAGS := --debug -DDEBUG
RELEASEFLAGS := -Ofast -march=native -mtune=native

release: $(SRC_FILES)
	clang++ $(CXXFLAGS) $(RELEASEFLAGS) $(LDFLAGS) -o $@ $^

debug: $(SRC_FILES)
	clang++ $(CXXFLAGS) $(DEBUGFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -rf release debug
