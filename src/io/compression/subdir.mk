# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
src/io/compression/bgzf_controller.cpp \
src/io/compression/tgzf_controller.cpp \
src/io/compression/tgzf_controller_stream.cpp 

OBJS += \
src/io/compression/bgzf_controller.o \
src/io/compression/tgzf_controller.o \
src/io/compression/tgzf_controller_stream.o 

CPP_DEPS += \
src/io/compression/bgzf_controller.d \
src/io/compression/tgzf_controller.d \
src/io/compression/tgzf_controller_stream.d 


# Each subdirectory must supply rules for building sources it contributes
src/io/compression/%.o: src/io/compression/%.cpp
	g++ -std=c++0x -I"src/" -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"

