# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/romanmarakulin/Projects/jpeg_decoder/jpeg_decoder

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/romanmarakulin/Projects/jpeg_decoder/jpeg_decoder/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/test_baseline.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/test_baseline.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test_baseline.dir/flags.make

CMakeFiles/test_baseline.dir/src/test_baseline.cpp.o: CMakeFiles/test_baseline.dir/flags.make
CMakeFiles/test_baseline.dir/src/test_baseline.cpp.o: ../src/test_baseline.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/romanmarakulin/Projects/jpeg_decoder/jpeg_decoder/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/test_baseline.dir/src/test_baseline.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_baseline.dir/src/test_baseline.cpp.o -c /Users/romanmarakulin/Projects/jpeg_decoder/jpeg_decoder/src/test_baseline.cpp

CMakeFiles/test_baseline.dir/src/test_baseline.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_baseline.dir/src/test_baseline.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/romanmarakulin/Projects/jpeg_decoder/jpeg_decoder/src/test_baseline.cpp > CMakeFiles/test_baseline.dir/src/test_baseline.cpp.i

CMakeFiles/test_baseline.dir/src/test_baseline.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_baseline.dir/src/test_baseline.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/romanmarakulin/Projects/jpeg_decoder/jpeg_decoder/src/test_baseline.cpp -o CMakeFiles/test_baseline.dir/src/test_baseline.cpp.s

CMakeFiles/test_baseline.dir/src/test_baseline.cpp.o.requires:

.PHONY : CMakeFiles/test_baseline.dir/src/test_baseline.cpp.o.requires

CMakeFiles/test_baseline.dir/src/test_baseline.cpp.o.provides: CMakeFiles/test_baseline.dir/src/test_baseline.cpp.o.requires
	$(MAKE) -f CMakeFiles/test_baseline.dir/build.make CMakeFiles/test_baseline.dir/src/test_baseline.cpp.o.provides.build
.PHONY : CMakeFiles/test_baseline.dir/src/test_baseline.cpp.o.provides

CMakeFiles/test_baseline.dir/src/test_baseline.cpp.o.provides.build: CMakeFiles/test_baseline.dir/src/test_baseline.cpp.o


CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.o: CMakeFiles/test_baseline.dir/flags.make
CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.o: /Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/romanmarakulin/Projects/jpeg_decoder/jpeg_decoder/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.o -c /Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp

CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp > CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.i

CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp -o CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.s

CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.o.requires:

.PHONY : CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.o.requires

CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.o.provides: CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.o.requires
	$(MAKE) -f CMakeFiles/test_baseline.dir/build.make CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.o.provides.build
.PHONY : CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.o.provides

CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.o.provides.build: CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.o


# Object files for target test_baseline
test_baseline_OBJECTS = \
"CMakeFiles/test_baseline.dir/src/test_baseline.cpp.o" \
"CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.o"

# External object files for target test_baseline
test_baseline_EXTERNAL_OBJECTS =

test_baseline: CMakeFiles/test_baseline.dir/src/test_baseline.cpp.o
test_baseline: CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.o
test_baseline: CMakeFiles/test_baseline.dir/build.make
test_baseline: libdecoder-lib.a
test_baseline: /usr/local/lib/libfftw3.dylib
test_baseline: /usr/local/lib/libpng.dylib
test_baseline: /usr/local/lib/libjpeg.dylib
test_baseline: CMakeFiles/test_baseline.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/romanmarakulin/Projects/jpeg_decoder/jpeg_decoder/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable test_baseline"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_baseline.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test_baseline.dir/build: test_baseline

.PHONY : CMakeFiles/test_baseline.dir/build

CMakeFiles/test_baseline.dir/requires: CMakeFiles/test_baseline.dir/src/test_baseline.cpp.o.requires
CMakeFiles/test_baseline.dir/requires: CMakeFiles/test_baseline.dir/Users/romanmarakulin/Projects/jpeg_decoder/contrib/catch_main.cpp.o.requires

.PHONY : CMakeFiles/test_baseline.dir/requires

CMakeFiles/test_baseline.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test_baseline.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test_baseline.dir/clean

CMakeFiles/test_baseline.dir/depend:
	cd /Users/romanmarakulin/Projects/jpeg_decoder/jpeg_decoder/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/romanmarakulin/Projects/jpeg_decoder/jpeg_decoder /Users/romanmarakulin/Projects/jpeg_decoder/jpeg_decoder /Users/romanmarakulin/Projects/jpeg_decoder/jpeg_decoder/cmake-build-debug /Users/romanmarakulin/Projects/jpeg_decoder/jpeg_decoder/cmake-build-debug /Users/romanmarakulin/Projects/jpeg_decoder/jpeg_decoder/cmake-build-debug/CMakeFiles/test_baseline.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test_baseline.dir/depend

