TARGET := router
SRCDIR := ./src
objDIR := ./obj

CXX ?= g++
CXXFLAGS := -std=c++17 -I${SRCDIR} -MMD -MP -Wall -Wextra -Wshadow -Wconversion
LDFLAGS  := -lpthread
RC_FLAGS := -Ofast
PF_FLAGS ?= -Ofast -pg -g
DC_FLAGS := -g -DDEBUG -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_GLIBCXX_DEBUG_PEDANTIC_ASSERT

ifeq ($(DEBUG), 1)
	CXXFLAGS += $(DC_FLAGS)
	OBJDIR   := $(objDIR)/debug
	TARGET   := $(TARGET).debug
else ifeq ($(PERF), 1)
	CXXFLAGS += $(PF_FLAGS)
	LDFLAGS  += $(PF_FLAGS)
	OBJDIR   := $(objDIR)/perf
	TARGET   := $(TARGET).perf
else
	CXXFLAGS += $(RC_FLAGS)
	OBJDIR   := $(objDIR)/release
endif

SRCS := $(shell find ${SRCDIR} -name '*.cpp')
OBJS := $(addprefix $(OBJDIR)/, $(patsubst ${SRCDIR}/%.cpp, %.o, $(SRCS)))
DEPS := $(addprefix $(OBJDIR)/, $(patsubst ${SRCDIR}/%.cpp, %.d, $(SRCS)))

IGNORE_OBJS := LayerAssignment.o
$(OBJDIR)/$(IGNORE_OBJS): CXXFLAGS += -w

$(shell [ -d $(OBJDIR) ] || mkdir -p $(OBJDIR))

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

-include $(DEPS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	$(RM) -r $(objDIR) $(TARGET){,.debug}
