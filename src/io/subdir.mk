# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
src/io/basic_writers.cpp \
src/io/output_writer.cpp \
src/io/reader.cpp 

OBJS += \
src/io/basic_writers.o \
src/io/output_writer.o \
src/io/reader.o 

CPP_DEPS += \
src/io/basic_writers.d \
src/io/output_writer.d \
src/io/reader.d 


# Each subdirectory must supply rules for building sources it contributes
src/io/%.o: src/io/%.cpp
	g++ $(CXXFLAGS) $(INCLUDE_PATH) -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"


