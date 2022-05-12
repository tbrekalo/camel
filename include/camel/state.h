#ifndef CAMEL_STATE_H_
#define CAMEL_STATE_H_

#include <filesystem>

#include "thread_pool/thread_pool.hpp"

namespace camel {

struct State {
  std::shared_ptr<thread_pool::ThreadPool> thread_pool;
  std::filesystem::path log_path;
};

};  // namespace camel

#endif /* CAMEL_STATE_H_ */
