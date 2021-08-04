#include "MinHeap.hpp"

void MinHeap::minHeapify(int idx)
{
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;

    if (left < size && array[left].key < array[smallest].key)
        smallest = left;

    if (right < size && array[right].key < array[smallest].key)
        smallest = right;

    if (smallest != idx)
    {
        // The nodes to be swapped in min heap
        MinHeapNode smallestNode = array[smallest];
        MinHeapNode idxNode = array[idx];

        // Swap positions
        pos[smallestNode.v] = idx;
        pos[idxNode.v] = smallest;

        MinHeapNode temp = array[smallest];
        array[smallest] = array[idx];
        array[idx] = temp;

        minHeapify(smallest);
    }
}

MinHeapNode MinHeap::extractMin()
{

    // Store the root node
    MinHeapNode root = array[0];

    // Replace root node with last node
    MinHeapNode lastNode = array[size - 1];
    array[0] = lastNode;

    // Update position of last node
    pos[root.v] = size - 1;
    pos[lastNode.v] = 0;

    // Reduce heap size and heapify root
    size--;
    minHeapify(0);

    return root;
}

void MinHeap::decreaseKey(int v, int key)
{
    // Get the index of v in  heap array
    int i = pos[v];

    // Get the node and update its key value
    array[i].key = key;

    // Travel up while the complete tree is not hepified.
    // This is a O(Logn) loop
    while (i && array[i].key < array[(i - 1) / 2].key)
    {
        // Swap this node with its parent
        pos[array[i].v] = (i - 1) / 2;
        pos[array[(i - 1) / 2].v] = i;

        MinHeapNode temp = array[i];
        array[i] = array[(i - 1) / 2];
        array[(i - 1) / 2] = temp;

        // move to parent index
        i = (i - 1) / 2;
    }
}
