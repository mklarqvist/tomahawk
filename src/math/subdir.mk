# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
src/math/fisher_math.cpp 

OBJS += \
src/math/fisher_math.o 

CPP_DEPS += \
src/math/fisher_math.d 


# Each subdirectory must supply rules for building sources it contributes
src/math/%.o: src/math/%.cpp
	g++ -std=c++0x -I"src/" -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"


