# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
src/main.cpp 

OBJS += \
src/main.o 

CPP_DEPS += \
src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: src/%.cpp
	g++ $(CXXFLAGS) $(INCLUDE_PATH) -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
