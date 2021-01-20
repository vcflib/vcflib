#include <memory>

namespace Detail {
  template <typename T, typename... Args>
  std::unique_ptr<T> makeUnique(Args &&... args) {
    return unique_ptr<T>(new T(forward<Args>(args)...));
  }
}
