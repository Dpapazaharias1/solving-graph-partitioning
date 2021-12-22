#ifndef _MinHeap_h
#define _MinHeap_h

#include <vector>

struct MinHeapNode
{
    int v;
    double key;
    inline MinHeapNode(int v, double key)
    {
        this->v = v;
        this->key = key;
    }
    inline MinHeapNode() {}
};

class MinHeap
{

public:
    int size;     /**< Number of heap nodes present currently */
    int capacity; /**< Capacity of min heap*/
    // This is needed for decreaseKey()
    std::vector<int> pos;
    std::vector<MinHeapNode> array;

    inline MinHeap(
        int capacity)
    {
        pos = std::vector<int>(capacity, 0);
        size = 0;
        this->capacity = capacity;
        array = std::vector<MinHeapNode>(capacity);
    };

    // A standard function to heapify at given idx
    // This function also updates position of nodes when they are swapped.
    // Position is needed for decreaseKey()
    void minHeapify(int idx);
    // Function to check if the given minHeap is empty or not
    inline int isEmpty() { return size == 0; };
    // Function to extract minimum node from heap
    MinHeapNode extractMin();

    // Function to decreasy key value of a given vertex v. This function
    // uses pos[] of min heap to get the current index of node in min heap
    void decreaseKey(int v, int key);

    // Function to check if a given vertex v is in min heap or not
    inline bool isInMinHeap(int v)
    {
        if (pos[v] < size)
            return true;
        return false;
    }
};
#endif
