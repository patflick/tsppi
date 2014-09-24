/*
 * Copyright (c) 2014 Patrick Flick <patrick.flick@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#ifndef TSPPI_PPI_UTILS_H
#define TSPPI_PPI_UTILS_H

#include <string>

/// The PPis used in this work
std::string ppis[5] = {"ccsb", "string", "bossi", "psicquic_all", "havu"};
/// The expression datasets used
std::string expr[5] = {"hpa", "hpa_all", "emtab", "rnaseq_atlas", "gene_atlas"};


/**
 * @brief Executes the given function for each PPI.
 */
template<typename T>
void foreach_ppi(T handle)
{
    for (std::string p : ppis)
    {
        handle(p);
    }
}

/**
 * @brief Executes the given function for all combinations of
 *        PPI networks and expression datasets.
 */
template<typename T>
void foreach_ppi_expr(T handle)
{
    for (std::string p : ppis)
        for (std::string e : expr)
        {
            handle(p, e);
        }
}

#endif // TSPPI_PPI_UTILS_H
