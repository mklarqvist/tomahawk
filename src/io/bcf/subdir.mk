# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
src/io/bcf/bcf_entry.cpp \
src/io/bcf/bcf_reader.cpp 

OBJS += \
src/io/bcf/bcf_entry.o \
src/io/bcf/bcf_reader.o 

CPP_DEPS += \
src/io/bcf/bcf_entry.d \
src/io/bcf/bcf_reader.d 


# Each subdirectory must supply rules for building sources it contributes
src/io/bcf/%.o: src/io/bcf/%.cpp
	g++ $(CXXFLAGS) $(INCLUDE_PATH) -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"


