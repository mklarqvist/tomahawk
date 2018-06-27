# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
src/tomahawk/tomahawk_calc.cpp \
src/tomahawk/tomahawk_importer.cpp \
src/tomahawk/tomahawk_reader.cpp \
src/tomahawk/import_writer.cpp \
src/tomahawk/tomahawk_output_reader.cpp \
src/tomahawk/output_entry.cpp \
src/tomahawk/output_entry_support.cpp \
src/tomahawk/output_filter.cpp 

OBJS += \
src/tomahawk/tomahawk_calc.o \
src/tomahawk/tomahawk_importer.o \
src/tomahawk/tomahawk_reader.o \
src/tomahawk/import_writer.o \
src/tomahawk/tomahawk_output_reader.o \
src/tomahawk/output_entry.o \
src/tomahawk/output_entry_support.o \
src/tomahawk/output_filter.o 

CPP_DEPS += \
src/tomahawk/tomahawk_calc.d \
src/tomahawk/tomahawk_importer.d \
src/tomahawk/tomahawk_reader.d \
src/tomahawk/import_writer.d \
src/tomahawk/tomahawk_output_reader.d \
src/tomahawk/output_entry.d \
src/tomahawk/output_entry_support.d \
src/tomahawk/output_filter.d 

# Each subdirectory must supply rules for building sources it contributes
src/tomahawk/%.o: src/tomahawk/%.cpp
	g++ $(CXXFLAGS) $(INCLUDE_PATH) -c -DVERSION=\"$(GIT_VERSION)\" -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"


