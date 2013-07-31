#ifndef Xfd92b845578be6e1eebb936ca29552e1ea3a2b24 
#define Xfd92b845578be6e1eebb936ca29552e1ea3a2b24 

#include <cassert>

template <typename CLASS>
class AlignedAllocator {
public:
    AlignedAllocator () {
        master_chunk.begin = master_chunk.end = 0;
        master_chunk.data = 0;
        master_chunk.next = 0;
    }

    CLASS *allocate () {
        return master_chunk.allocate ();
    }

    static CLASS *find_object (void *pointer) {
        return Chunk::find_element (pointer)->payload ();
    }

    static const CLASS *find_object (const void *pointer) {
        return const_cast <const CLASS *> (
            find_object (const_cast <void *> (pointer)));
    }

private:
    struct Chunk {
        ~Chunk () {
            clear ();
        }
        
        void clear () {
            for (Element *pos = find_begin_pointer (data); pos != begin; ++pos)
                pos->~Element ();
            begin = end = 0;
            delete[] data;
            data = 0;
            delete next;
            next = 0;
        }

        class Element {
        public:
            CLASS stuff;
            Element () {
                assert (misalignment (this) == 0);
#ifndef NDEBUG
                magic = this;
#endif
            }

            ~Element () {
                assert (misalignment (this) == 0);
                assert (magic == this);
#ifndef NDEBUG
                magic = 0;
#endif
            }

            inline CLASS *payload () {
                assert (this->magic == this);
                return &stuff;
            }

        private:
#ifndef NDEBUG
            Element *magic;
#endif
            Element (const Element &);
            void operator= (const Element &);
        };

        static size_t sizeof_element () {
            return reinterpret_cast <char *> (static_cast <Element *> (0) + 1)
                - reinterpret_cast <char *> (0);
        }

        static size_t misalignment (void *pointer) {
            size_t addr = reinterpret_cast <size_t> (pointer);
            return addr % sizeof_element ();
        }

        static Element *find_element (void *pointer) {
            char *begin = reinterpret_cast <char *> (pointer)
                - misalignment (pointer);
            Element *ret = reinterpret_cast <Element *> (begin);
            assert (reinterpret_cast <char *> (pointer) >= reinterpret_cast <char *> (ret));
            return ret;
        }

        static Element *find_begin_pointer (char *memblock_) {
            char *memblock = memblock_ + sizeof_element () - 1;
            Element *ret = find_element (memblock);
            assert (reinterpret_cast <char *> (ret) >= memblock_);
            return ret;
        }

        CLASS *allocate () {
            if (begin != end) {
                new (begin) Element ();
                return (begin++)->payload ();
            } else {
                // we need to allocate a new Chunk first
                this->next = new Chunk (*this);
                this->data = 0;
                this->begin = 0;
                this->end = 0;
                this->data = new char [100000*sizeof_element ()];
                // align properly
                begin = find_begin_pointer (this->data);
                end = begin + 100000 - 1;
                // now try again
                return allocate ();
            }
        }

        char *data;
        Element *begin, *end;
        Chunk *next;
    }
    master_chunk;
};

#endif // Xfd92b845578be6e1eebb936ca29552e1ea3a2b24
