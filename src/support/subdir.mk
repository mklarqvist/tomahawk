# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
src/support/magic_constants.cpp \
src/support/helpers.cpp 

OBJS += \
src/support/magic_constants.o \
src/support/helpers.o 

CPP_DEPS += \
src/support/magic_constants.d \
src/support/helpers.d 


# Each subdirectory must supply rules for building sources it contributes
src/support/%.o: src/support/%.cpp
	g++ $(CXXFLAGS) $(INCLUDE_PATH) -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"


