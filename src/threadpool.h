#ifndef _THREAD_POOL_H_
#define _THREAD_POOL_H_

#include <atomic>
#include <condition_variable>
#include <functional>
#include <future>
#include <memory>
#include <map>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

template<typename T>
class threadsafe_queue {
private:
	mutable std::mutex mut;
	std::queue<T> data_queue;

public:
	threadsafe_queue() {}

	void push(T &&new_value) {
		std::lock_guard<std::mutex> lk(mut);
		data_queue.emplace(new_value);
	}

	bool try_pop(T& value) {
		std::lock_guard<std::mutex> lk(mut);
		if (data_queue.empty())
			return false;
		value = std::move(data_queue.front());
		data_queue.pop();
		return true;
	}

	size_t size() {
		std::lock_guard<std::mutex> lk(mut);
		return data_queue.size();
	}
};

class thread_pool
{
	std::atomic_bool done;
	int nthreads;
        std::map<std::thread::id, int> thread_id;
	threadsafe_queue<std::function<void()>> work_queue;
	std::vector<std::thread> threads;
        std::condition_variable cv;
        std::mutex m;

	void worker_thread() {
		while (!done) {
			std::function<void()> task;
                        if (work_queue.try_pop(task)) {
				task();
                        } else {
                                std::unique_lock<std::mutex> lock(m);
                                cv.wait(lock, [&] {return work_queue.size() != 0 || done; });
                        }
                }
	}
public:
	thread_pool(int nthr):
		done(false), nthreads(nthr)
		{
			try {
                                for (int i = 0; i < nthreads; ++i) {
                                        threads.emplace_back(
                                                std::thread(&thread_pool::worker_thread, this));
                                        thread_id[threads[i].get_id()] = i;
                                }
                        } catch (...) {
				done = true;
				throw;
			}
		}
	~thread_pool() {
                while (work_queue.size() != 0) ;
		done = true;
                std::unique_lock<std::mutex> lock(m);
                cv.notify_all();
                lock.unlock();
		for (std::thread& thread: threads) {
			thread.join();
		}

	}

	template<typename Function, typename... Args>
	std::future<typename std::result_of<Function(Args...)>::type>
	submit(Function &&f, Args&&... args) {
		using result_type = typename std::result_of<Function(Args...)>::type;
		auto task = std::make_shared<std::packaged_task<result_type()>>(std::bind(std::forward<Function>(f), std::forward<Args>(args)...));
		std::future<result_type> res(task->get_future());
                work_queue.push([task] { (*task)(); });
                std::unique_lock<std::mutex> lock(m);
                cv.notify_one();
                return res;
	}

        int size() {
                return nthreads;
        }

        int thread_id_to_int(std::thread::id id) {
                return thread_id[id];
        }

        template<typename T, typename Function>
        void parallel_for(T start, T end, T stride, Function &&f) {
                T range = end - start;
                T block_size = range / (nthreads);
                T block_start = start;
                T block_end = block_start + block_size;
                if (block_size == 0)
                        block_end = end;
                std::vector<std::future<void>> res;
                while (block_start < end) {
                        res.emplace_back(submit(f, block_start, block_end, stride));
                        block_start = block_end;
                        block_end = block_end + block_size;
                        if (block_end >= end)
                                block_end = end;
                }
                for (size_t i = 0; i < res.size(); i++)
                        res[i].get();
        }
};

#endif
