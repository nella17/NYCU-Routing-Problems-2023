TARGET := router
SRCDIR := ./src
objDIR := ./obj

CXX ?= g++
CXXFLAGS := -std=c++17 -I${SRCDIR} -MMD -MP
LDFLAGS  := 
RC_FLAGS := -Ofast
DC_FLAGS := -g -DDEBUG -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_GLIBCXX_DEBUG_PEDANTIC_ASSERT

ifeq ($(DEBUG), 1)
	CXXFLAGS += $(DC_FLAGS)
	OBJDIR   := $(objDIR)/debug
	TARGET   := $(TARGET).debug
else
	CXXFLAGS += $(RC_FLAGS)
	OBJDIR   := $(objDIR)/release
endif

SRCS := $(shell find ${SRCDIR} -name '*.cpp')
OBJS := $(addprefix $(OBJDIR)/, $(patsubst ${SRCDIR}/%.cpp, %.o, $(SRCS)))
DEPS := $(addprefix $(OBJDIR)/, $(patsubst ${SRCDIR}/%.cpp, %.d, $(SRCS)))

$(shell [ -d $(OBJDIR) ] || mkdir -p $(OBJDIR))

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

-include $(DEPS)

$(OBJDIR)%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	$(RM) -r $(objDIR) $(TARGET){,.debug}
