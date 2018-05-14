# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
src/io/vcf/vcf_header.cpp 

OBJS += \
src/io/vcf/vcf_header.o 

CPP_DEPS += \
src/io/vcf/vcf_header.d 


# Each subdirectory must supply rules for building sources it contributes
src/io/vcf/%.o: src/io/vcf/%.cpp
	g++ -std=c++0x -I"src/" -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"


