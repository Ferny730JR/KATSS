find_package(Threads REQUIRED)

include(CheckCSourceCompiles)

set(CMAKE_REQUIRED_LIBRARIES Threads::Threads)
check_c_source_compiles("
#include <threads.h>

int main(void)
{
	mtx_t mutex;

	if (mtx_init(&mutex, mtx_plain) != thrd_success)
		return 1;
	if (mtx_lock(&mutex) != thrd_success)
		return 2;
	if (mtx_unlock(&mutex) != thrd_success)
		return 3;
	mtx_destroy(&mutex);
	return 0;
}
" SEQF_HAVE_C11_THREADS)
unset(CMAKE_REQUIRED_LIBRARIES)

add_library(seqf_thread_config INTERFACE)
target_compile_definitions(seqf_thread_config INTERFACE
	SEQF_USE_C11_THREADS=$<BOOL:${SEQF_HAVE_C11_THREADS}>)
target_link_libraries(seqf_thread_config INTERFACE Threads::Threads)

if(NOT SEQF_HAVE_C11_THREADS)
	add_subdirectory(tinycthread)
endif()
