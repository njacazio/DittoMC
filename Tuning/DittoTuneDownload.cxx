///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   DittoTuneDownload.cxx
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/19
/// \brief  Download and cache remote Ditto tune files.
///

#include "DittoTuneDownload.h"

#include <curl/curl.h>

#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/file.h>
#include <sys/types.h>

#include <fcntl.h>
#include <unistd.h>
#endif

namespace Ditto::Tunes
{
namespace
{

class CurlGlobal
{
 public:
  CurlGlobal()
  {
    const CURLcode code = curl_global_init(CURL_GLOBAL_DEFAULT);
    if (code != CURLE_OK) {
      throw std::runtime_error(std::string("Ditto: curl_global_init failed: ") +
                               curl_easy_strerror(code));
    }
  }

  ~CurlGlobal() { curl_global_cleanup(); }
};

void ensureCurlInitialized()
{
  static const CurlGlobal global;
  (void)global;
}

std::filesystem::path temporaryRoot()
{
  constexpr std::array<const char*, 3> variables = {"TMP", "TMPDIR", "TEMP"};
  for (const char* variable : variables) {
    if (const char* value = std::getenv(variable); value && *value) {
      return value;
    }
  }
  return std::filesystem::temp_directory_path();
}

std::uint64_t fnv1a64(std::string_view text)
{
  std::uint64_t hash = 14695981039346656037ULL;
  for (const unsigned char c : text) {
    hash ^= static_cast<std::uint64_t>(c);
    hash *= 1099511628211ULL;
  }
  return hash;
}

std::string hexadecimal(std::uint64_t value)
{
  std::ostringstream stream;
  stream << std::hex << std::setfill('0') << std::setw(16) << value;
  return stream.str();
}

std::string fileNameFromUrl(std::string_view url)
{
  const std::size_t end = url.find_first_of("?#");
  const std::string_view path = url.substr(0, end);
  const std::size_t slash = path.find_last_of('/');
  std::string name(path.substr(slash == std::string_view::npos ? 0 : slash + 1));

  if (name.empty()) {
    throw std::invalid_argument("Ditto: tune URL has no file name: " + std::string(url));
  }

  for (char& c : name) {
    const bool allowed =
      (c >= 'a' && c <= 'z') ||
      (c >= 'A' && c <= 'Z') ||
      (c >= '0' && c <= '9') || c == '.' || c == '_' || c == '-';
    if (!allowed) {
      c = '_';
    }
  }
  return name;
}

size_t writeToFile(char* data, size_t size, size_t count, void* userData)
{
  return std::fwrite(data, size, count, static_cast<std::FILE*>(userData));
}

std::filesystem::path makePartPath(const std::filesystem::path& destination)
{
  std::string suffix = ".part.";
#if defined(__unix__) || defined(__APPLE__)
  suffix += std::to_string(static_cast<long long>(::getpid()));
#else
  suffix += std::to_string(
    std::chrono::steady_clock::now().time_since_epoch().count());
#endif
  return destination.string() + suffix;
}

class DownloadLock
{
 public:
  explicit DownloadLock(const std::filesystem::path& path)
  {
#if defined(__unix__) || defined(__APPLE__)
    mFd = ::open(path.c_str(), O_CREAT | O_RDWR, 0600);
    if (mFd < 0) {
      throw std::runtime_error("Ditto: cannot create tune-cache lock " + path.string());
    }
    if (::flock(mFd, LOCK_EX) != 0) {
      ::close(mFd);
      mFd = -1;
      throw std::runtime_error("Ditto: cannot lock tune cache " + path.string());
    }
#else
    (void)path;
#endif
  }

  DownloadLock(const DownloadLock&) = delete;
  DownloadLock& operator=(const DownloadLock&) = delete;

  ~DownloadLock()
  {
#if defined(__unix__) || defined(__APPLE__)
    if (mFd >= 0) {
      ::flock(mFd, LOCK_UN);
      ::close(mFd);
    }
#endif
  }

 private:
#if defined(__unix__) || defined(__APPLE__)
  int mFd = -1;
#endif
};

void performDownload(std::string_view url, const std::filesystem::path& destination)
{
  ensureCurlInitialized();

  const std::filesystem::path partPath = makePartPath(destination);
  std::error_code ec;
  std::filesystem::remove(partPath, ec);

  struct FileCloser {
    void operator()(std::FILE* file) const noexcept
    {
      if (file) {
        std::fclose(file);
      }
    }
  };
  std::unique_ptr<std::FILE, FileCloser> output(std::fopen(partPath.c_str(), "wb"));
  if (!output) {
    throw std::runtime_error("Ditto: cannot create temporary tune file " +
                             partPath.string());
  }

  CURL* rawCurl = curl_easy_init();
  if (!rawCurl) {
    std::filesystem::remove(partPath, ec);
    throw std::runtime_error("Ditto: curl_easy_init failed");
  }
  std::unique_ptr<CURL, decltype(&curl_easy_cleanup)> curl(rawCurl, &curl_easy_cleanup);

  const std::string urlString(url);
  std::array<char, CURL_ERROR_SIZE> errorBuffer{};

  curl_easy_setopt(curl.get(), CURLOPT_URL, urlString.c_str());
  curl_easy_setopt(curl.get(), CURLOPT_FOLLOWLOCATION, 1L);
  curl_easy_setopt(curl.get(), CURLOPT_FAILONERROR, 1L);
  curl_easy_setopt(curl.get(), CURLOPT_NOSIGNAL, 1L);
  curl_easy_setopt(curl.get(), CURLOPT_CONNECTTIMEOUT, 30L);
  curl_easy_setopt(curl.get(), CURLOPT_LOW_SPEED_LIMIT, 1L);
  curl_easy_setopt(curl.get(), CURLOPT_LOW_SPEED_TIME, 60L);
  curl_easy_setopt(curl.get(), CURLOPT_USERAGENT, "DittoMC tune downloader");
  curl_easy_setopt(curl.get(), CURLOPT_ERRORBUFFER, errorBuffer.data());
  curl_easy_setopt(curl.get(), CURLOPT_WRITEFUNCTION, &writeToFile);
  curl_easy_setopt(curl.get(), CURLOPT_WRITEDATA, output.get());

  const CURLcode result = curl_easy_perform(curl.get());
  output.reset();

  if (result != CURLE_OK) {
    std::filesystem::remove(partPath, ec);
    std::string message = "Ditto: failed to download tune from " + urlString + ": ";
    message += errorBuffer[0] ? errorBuffer.data() : curl_easy_strerror(result);
    throw std::runtime_error(message);
  }

  const auto size = std::filesystem::file_size(partPath, ec);
  if (ec || size == 0) {
    std::filesystem::remove(partPath, ec);
    throw std::runtime_error("Ditto: downloaded tune is empty: " + urlString);
  }

  // The cache lock makes the destination private to this downloader while the
  // rename happens. rename() is atomic on the local filesystems normally used
  // for TMPDIR/TMP, so readers never observe an incomplete ROOT file.
  std::filesystem::rename(partPath, destination, ec);
  if (ec) {
    std::filesystem::remove(partPath, ec);
    throw std::runtime_error("Ditto: cannot move downloaded tune into cache: " +
                             destination.string());
  }
}

} // namespace

std::filesystem::path cacheDirectory()
{
  return temporaryRoot() / "DittoCache" / "tunes";
}

bool isRemote(std::string_view source)
{
  return source.rfind("https://", 0) == 0 || source.rfind("http://", 0) == 0;
}

std::filesystem::path cachedPath(std::string_view url)
{
  if (!isRemote(url)) {
    throw std::invalid_argument("Ditto: cachedPath requires an HTTP(S) URL");
  }

  // Include a URL hash in the path so identically named assets from different
  // releases/servers cannot alias each other in the cache.
  return cacheDirectory() / hexadecimal(fnv1a64(url)) / fileNameFromUrl(url);
}

std::filesystem::path fetch(std::string_view url, bool force)
{
  if (!isRemote(url)) {
    throw std::invalid_argument("Ditto: fetch requires an HTTP(S) URL: " + std::string(url));
  }

  const std::filesystem::path destination = cachedPath(url);
  std::error_code ec;
  std::filesystem::create_directories(destination.parent_path(), ec);
  if (ec) {
    throw std::runtime_error("Ditto: cannot create tune cache directory " + destination.parent_path().string() + ": " + ec.message());
  }

  // Serialize downloads of the same URL across processes sharing the job TMP.
  DownloadLock lock(destination.string() + ".lock");

  if (!force && std::filesystem::is_regular_file(destination, ec) && !ec &&
      std::filesystem::file_size(destination, ec) > 0 && !ec) {
    return destination;
  }

  // Remove a forced or invalid cache entry before the atomic rename.
  if (force || std::filesystem::exists(destination, ec)) {
    std::filesystem::remove(destination, ec);
    ec.clear();
  }

  std::cerr << "Ditto: downloading tune " << url << " -> " << destination << '\n';
  performDownload(url, destination);
  return destination;
}

std::filesystem::path resolve(std::string_view source, bool force)
{
  if (source.empty()) {
    throw std::invalid_argument("Ditto: empty tune location");
  }

  if (isRemote(source)) {
    return fetch(source, force);
  }

  const std::filesystem::path path(source);
  std::error_code ec;
  if (!std::filesystem::is_regular_file(path, ec) || ec) {
    throw std::runtime_error("Ditto: tune file does not exist: " + path.string());
  }
  return path;
}

} // namespace Ditto::Tunes
