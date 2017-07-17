#ifndef SPINLOCK_H_
#define SPINLOCK_H_

#include <atomic>

namespace Tomahawk{
namespace Algorithm{

class SpinLock{
public:
    void lock(){
        while(lck.test_and_set(std::memory_order_acquire))
        {}
    }
    void unlock(){lck.clear(std::memory_order_release);}

private:
    std::atomic_flag lck = ATOMIC_FLAG_INIT;
};

}
}

#endif /* SPINLOCK_H_ */
