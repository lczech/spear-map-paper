/*
    SPEAR - Sorting Petabytes of Environmental and Ancient Reads
    Copyright (C) 2026 Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@sund.ku.dk>
    University of Copenhagen, Globe Institute, Section for GeoGenetics
    Oster Voldgade 5-7, 1350 Copenhagen K, Denmark
*/

#pragma once

#include <coroutine>
#include <exception>
#include <utility>

// =================================================================================================
//     Generator
// =================================================================================================

// A minimal C++20 coroutine generator that allows lazy, streaming iteration via range-based for.
// Usage:
//     Generator<T> producer() { co_yield value; }
//     for( auto const& v : producer() ) { ... }
template<typename T>
class Generator
{
public:

    // -------------------------------------------------------------------------
    //     Promise
    // -------------------------------------------------------------------------

    struct promise_type
    {
        T value;

        Generator get_return_object() {
            return Generator{ std::coroutine_handle<promise_type>::from_promise( *this ) };
        }
        std::suspend_always initial_suspend() noexcept { return {}; }
        std::suspend_always final_suspend()   noexcept { return {}; }
        std::suspend_always yield_value( T v ) {
            value = std::move( v );
            return {};
        }
        void return_void() noexcept {}
        void unhandled_exception() { std::rethrow_exception( std::current_exception() ); }
    };

    // -------------------------------------------------------------------------
    //     Iterator
    // -------------------------------------------------------------------------

    struct Sentinel {};

    struct Iterator
    {
        std::coroutine_handle<promise_type> handle;

        Iterator& operator++() {
            handle.resume();
            return *this;
        }
        T const& operator*()  const noexcept { return handle.promise().value; }
        T const* operator->() const noexcept { return &handle.promise().value; }
        bool operator==( Sentinel ) const noexcept { return handle.done(); }
    };

    // -------------------------------------------------------------------------
    //     Generator
    // -------------------------------------------------------------------------

    explicit Generator( std::coroutine_handle<promise_type> h ) noexcept : handle_( h ) {}
    ~Generator() { if( handle_ ) handle_.destroy(); }

    Generator( Generator const& ) = delete;
    Generator& operator=( Generator const& ) = delete;

    Generator( Generator&& o ) noexcept : handle_( o.handle_ ) { o.handle_ = nullptr; }
    Generator& operator=( Generator&& o ) noexcept {
        if( this != &o ) {
            if( handle_ ) handle_.destroy();
            handle_ = o.handle_;
            o.handle_ = nullptr;
        }
        return *this;
    }

    // Resumes to the first co_yield before returning the iterator.
    Iterator begin() {
        handle_.resume();
        return Iterator{ handle_ };
    }
    Sentinel end() noexcept { return {}; }

private:
    std::coroutine_handle<promise_type> handle_ = nullptr;
};
