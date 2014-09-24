/*
 * Copyright (c) 2014 Patrick Flick <patrick.flick@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#ifndef TSPPI_UTILS_H
#define TSPPI_UTILS_H

// STL includes
#include <map>

/**
 * @brief Helper function to reverse a std::map and return it as a new map.
 *
 * This function reverses the roles of the Key and the Value of a map and creates
 * a new map. This should only be used for std::maps that are a 1:1 (one-to-one)
 * mapping. In case a value appears more than once, the behaviour is not specified.
 *
 * @param map   The map to be reversed.
 * @return      A new reversed map.
 */
template <typename S, typename T>
std::map<T, S> reverse_map(const std::map<S, T>& map)
{
    std::map<T, S> result;
    for (auto it : map)
    {
        result[it.second] = it.first;
    }
    return result;
}

#endif // #define TSPPI_UTILS_H
