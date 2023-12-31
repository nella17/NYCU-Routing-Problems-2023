CXX ?= g++
TARGET := Lab1
SRCDIR := .
OBJDIR := ./obj
CXXFLAGS := -std=c++17 -Wall -Wextra -Wshadow -Wconversion
RC_FLAGS := -Ofast
DC_FLAGS := -g -DDEBUG -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_GLIBCXX_DEBUG_PEDANTIC_ASSERT

ifeq ($(DEBUG), 1)
	CXXFLAGS := $(CXXFLAGS) $(DC_FLAGS)
	OBJDIR := $(OBJDIR)/debug
else
	CXXFLAGS := $(CXXFLAGS) $(RC_FLAGS)
	OBJDIR := $(OBJDIR)/release
endif

SRCS := $(notdir $(wildcard ${SRCDIR}/*.cpp))
OBJS := $(addprefix $(OBJDIR)/, $(patsubst %.cpp, %.o, $(SRCS)))
DEPS := $(addprefix $(OBJDIR)/, $(patsubst %.cpp, %.d, $(SRCS)))

$(shell [ -d $(OBJDIR) ] || mkdir -p $(OBJDIR))

all: $(TARGET) $(OBJS)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

-include $(DEPS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -MP -MMD -c $< -o $@

.PHONY: clean
clean:
	rm -rf $(OBJDIR) $(TARGET)
