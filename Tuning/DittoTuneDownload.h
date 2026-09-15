///
/// SPDX-License-Identifier: GPL-3.0-or-later
///
/// Copyright (C) 2026 Nicolò Jacazio
///
/// See LICENSE for details.
///
/// \file   DittoTuneDownload.h
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026/09/19
/// \brief  Download and cache remote Ditto tune files.
///

#pragma once

#include <filesystem>
#include <string_view>

namespace Ditto::Tunes
{

/// Return the job-local Ditto tune cache directory.
///
/// The base temporary directory is selected from TMP, TMPDIR, TEMP, and finally
/// std::filesystem::temp_directory_path(). The actual cache lives in
/// <tmp>/Ditto/tunes.
std::filesystem::path cacheDirectory();

/// Return true for HTTP(S) tune locations.
bool isRemote(std::string_view source);

/// Return the deterministic cache path associated with a remote URL.
/// No network operation is performed.
std::filesystem::path cachedPath(std::string_view url);

/// Download a remote tune into the temporary cache and return its local path.
///
/// If the tune is already cached, no network request is made unless force=true.
/// Downloads are written to a temporary file and atomically renamed into place,
/// so readers never observe a partial tune.
std::filesystem::path fetch(std::string_view url, bool force = false);

/// Resolve either a local tune path or an HTTP(S) URL to a local file.
///
/// Local paths are returned unchanged after checking that they exist. Remote
/// URLs are downloaded through fetch().
std::filesystem::path resolve(std::string_view source, bool force = false);

} // namespace Ditto::Tunes
