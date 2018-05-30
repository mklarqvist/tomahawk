# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
src/tomahawk/TomahawkCalc.cpp \
src/tomahawk/tomahawk_importer.cpp \
src/tomahawk/tomahawk_reader.cpp \
src/tomahawk/import_writer.cpp 

OBJS += \
src/tomahawk/TomahawkCalc.o \
src/tomahawk/tomahawk_importer.o \
src/tomahawk/tomahawk_reader.o \
src/tomahawk/import_writer.o 

CPP_DEPS += \
src/tomahawk/TomahawkCalc.d \
src/tomahawk/tomahawk_importer.d \
src/tomahawk/tomahawk_reader.d \
src/tomahawk/import_writer.d 

# Each subdirectory must supply rules for building sources it contributes
src/tomahawk/%.o: src/tomahawk/%.cpp
	g++ $(CXXFLAGS) $(INCLUDE_PATH) -c -DVERSION=\"$(GIT_VERSION)\" -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"


