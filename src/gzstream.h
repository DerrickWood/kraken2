#ifndef __GZSTREAM_H__
#define __GZSTREAM_H__

#include <assert.h>
#include <zlib.h>
#include <ios>
#include <locale>
#include <streambuf>
#include <iostream>
#include <string>
#include <vector>

template<typename charT, typename traits = std::char_traits<char>>
class basic_gzbuf : public std::basic_streambuf<charT, traits> {
        typedef charT char_type;
        typedef traits traits_type;
        typedef typename traits::int_type int_type;
        typedef typename traits::pos_type pos_type;
        typedef typename traits::off_type off_type;

public:
        basic_gzbuf(const char *filename);
        basic_gzbuf(const std::vector<std::string> *filenames);
        void open(const char *filename);
        void close();
        virtual ~basic_gzbuf() {
                close();
        }


protected:
        basic_gzbuf() {
                delete[] buffer_;
        }
        virtual int_type underflow();
        // virtual basic_gzbuf *setbuf(char_type *buf, std::streamsize size);

private:
        char_type *buffer_;
        std::streamsize size_;
        const std::vector<std::string> *filenames_;
        size_t file_idx_;
        static const int DEFAULT_BUFSIZ = 1024 * 8;
        gzFile f;
};

template<typename charT, typename traits>
basic_gzbuf<charT, traits>::basic_gzbuf(const char *filename)
        : buffer_(new char_type[DEFAULT_BUFSIZ]),
          size_(DEFAULT_BUFSIZ),
          filenames_(NULL)
{
        this->open(filename);
        this->setg(buffer_, buffer_ + size_, buffer_ + size_);
}

template<typename charT, typename traits>
basic_gzbuf<charT, traits>::basic_gzbuf(const std::vector<std::string> *filenames)
        : buffer_(new char_type[DEFAULT_BUFSIZ]),
          size_(DEFAULT_BUFSIZ),
          filenames_(filenames),
          file_idx_(0)
{
        assert(!filenames->empty());
        this->open(filenames_->at(file_idx_).c_str());
        file_idx_++;
        this->setg(buffer_, buffer_ + size_, buffer_ + size_);
}

template<typename charT, typename traits>
typename basic_gzbuf<charT, traits>::int_type
basic_gzbuf<charT, traits>::underflow() {
        int read = gzread(f, buffer_, size_);
        if (read == 0 && filenames_ != NULL && file_idx_ < filenames_->size()) {
                close();
                open(filenames_->at(file_idx_++).c_str());
                read = gzread(f, buffer_, size_);
        }
        this->setg(buffer_, buffer_, buffer_ + read);
        if (this->egptr() == this->gptr())
                return traits::eof();
        else
                return traits::to_int_type(*this->gptr());
}

template<typename charT, typename traits>
void basic_gzbuf<charT, traits>::open(const char *filename) {
        f = gzopen(filename, "rb");
        if (f == NULL) {
                std::cerr << "Unable to open \"" << filename << "\" for reading"
                          << std::endl;
                exit(1);
        }
        gzbuffer(f, DEFAULT_BUFSIZ);
}

template<typename charT, typename traits>
void basic_gzbuf<charT, traits>::close() {
        if (f != NULL)
                gzclose(f);
}

template<typename charT, typename traits = std::char_traits<charT>>
class basic_gzistream : public std::basic_istream<charT, traits>
{
public:
        typedef charT char_type;
        typedef traits traits_type;
        typedef typename traits::int_type int_type;
        typedef typename traits::pos_type pos_type;
        typedef typename traits::off_type off_type;

        basic_gzistream(const char *filename) : std::basic_istream<charT, traits>(NULL),
                                                sb(filename)
        {
                this->init(&sb);
        }

        basic_gzistream(const std::vector<std::string> *filenames): std::basic_istream<charT, traits>(NULL),
                                                              sb(filenames)
        {
                this->init(&sb);
        }
private:
        basic_gzbuf<char_type, traits_type> sb;
};

typedef basic_gzistream<char> gzistream;
#endif
