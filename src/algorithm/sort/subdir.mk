# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
src/algorithm/sort/output_sorter.cpp 

OBJS += \
src/algorithm/sort/output_sorter.o 

CPP_DEPS += \
src/algorithm/sort/output_sorter.d 


# Each subdirectory must supply rules for building sources it contributes
src/algorithm/sort/%.o: src/algorithm/sort/%.cpp
	g++ -std=c++0x -I"src/" -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"

