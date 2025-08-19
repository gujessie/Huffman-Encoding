#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <functional>
#include <queue>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <algorithm>   // for std::copy_n
#include "Floats_Huffman_Encoding.h"  // supplies MAXLENGTH & helpers
#include <gpgme.h>
#include "t-support.h"

using namespace std;

int get_frq_count(float* src, size_t num_elements, float error, float* min, float* max) {
    // Find min and max values in data
    *min = src[0];
    *max = src[0];
    for (size_t i = 1; i < num_elements; i++) {
        if (src[i] < *min)
            *min = src[i];
        if (src[i] > *max)
            *max = src[i];
    }

    // Calculate the interval and bucket width
    float interval = *max - *min;

    int frq_count = static_cast<int>(interval / (error * 2) + 1); //frq_count as in number of buckets
    float bucket_size = interval / frq_count;

    return frq_count;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input_file> <output_file>\n", argv[0]);
        return 1;
    }

    /* ---- open & load input ---- */
    FILE *input = fopen(argv[1], "rb");
    if (!input) {
        perror("Error opening input file");
        return 1;
    }

    fseek(input, 0, SEEK_END);
    const size_t num_elements = ftell(input) / sizeof(float);
    fseek(input, 0, SEEK_SET);

    printf("num_elements: %zu\n", num_elements);

    float *src = (float *)malloc(num_elements * sizeof(float));
    if (!src) {
        fprintf(stderr, "Error allocating memory for src\n");
        fclose(input);
        return 1;
    }

    const size_t read_elements = fread(src, sizeof(float), num_elements, input);
    fclose(input);
    if (read_elements != num_elements) {
        fprintf(stderr, "Error reading input file\n");
        free(src);
        return 1;
    }

    /* -----------------------------------------------------------

       Destination buffer sizing
       ----------------------------------------------------------- */
    const float error = 0.25f;
    float min, max;
    int frq_count = get_frq_count(src, num_elements, error, &min, &max);
    
    //Gary
    const size_t header = sizeof(size_t)          /* num_elements   */
    + 2 * sizeof(float)       /* min, bucket    */
    + sizeof(int)            /* frq_count      */
    +  frq_count * sizeof(int);     /* worst-case frqs[]              */

    const size_t dest_capacity =
      header
    + num_elements * sizeof(int);              /* worst-case 32-bit code each    */

    printf("dest_capacity = %lu, frq_count = %d\n", dest_capacity, frq_count);

    char *dest = (char *)malloc(dest_capacity);

    // ORIGINAL (fixed MAXLENGTH):
    // char *dest = (char *)malloc(MAXLENGTH);

    if (!dest) {
        perror("Error allocating memory for dest");
        free(src);
        return 1;
    }

    int destsize_bytes = 0;

    // after this call, the dest has ?
    if (compress(src, dest, num_elements, &destsize_bytes, error)) {
        fprintf(stderr, "Compression failed\n");
        free(src);
        free(dest);
        return 1;
    }

    //TODO: why dest_capacity is smaller destsize_bytes.
    printf("compressed + header: %d bytes (allocated %zu)\n", destsize_bytes, dest_capacity);

    // Debug: dump compressed+header bytes to stdout
    
    fwrite(dest, destsize_bytes, 1, stdout);

    // dest is compressed data. out is encrypted dest.
    // if we decrypt out, we should get dest back.

    //START ENCRYPTION HERE, DEST HAS COMPRESSED DATAf

 
    // begining of encryption
    gpgme_ctx_t ctx;
    gpgme_error_t err;
    gpgme_data_t in, out, decrpt_data;
    gpgme_key_t key[3] = { NULL, NULL, NULL };
    gpgme_encrypt_result_t result;

    init_gpgme (GPGME_PROTOCOL_OpenPGP);

    err = gpgme_new (&ctx);
    fail_if_err (err);
    gpgme_set_armor (ctx, 1);

    printf("JG E 2\n");

    /* Generate test keys. */
    char *fprs[2];
    err = generate_test_keys (ctx,2,fprs,NULL);
    fail_if_err (err);

    err = gpgme_data_new_from_mem (&in, dest, destsize_bytes, 0);
    fail_if_err (err);

    err = gpgme_data_new (&out);
    fail_if_err (err);

    err = gpgme_get_key (ctx,fprs[0],&key[0],0);
    fail_if_err (err);
    err = gpgme_get_key (ctx,fprs[1],&key[1],0);
    fail_if_err (err);

    //encryption here
    err = gpgme_op_encrypt (ctx, key, GPGME_ENCRYPT_ALWAYS_TRUST, in, out);
    fail_if_err (err);
    result = gpgme_op_encrypt_result (ctx);
    if (result->invalid_recipients)
    {
    fprintf (stderr, "Invalid recipient encountered: %s\n",
       result->invalid_recipients->fpr);
        exit (1);
    }

    //print encrypted data
    print_data (out);
    //printf("JG 3 E\n");

#if 1
    // decrypt

    gpgme_data_release (in);
    err = gpgme_data_new (&in);
    fail_if_err (err);

    /* Decrypt data. */
    // potentially the bug.
    err = gpgme_op_decrypt (ctx,out,in);
    fail_if_err (err);
    gpgme_decrypt_result_t dec_result;
    dec_result = gpgme_op_decrypt_result (ctx);
    if (dec_result->unsupported_algorithm) {
    fprintf (stderr, "%s:%i: unsupported algorithm: %s\n",
           __FILE__, __LINE__, dec_result->unsupported_algorithm);
    exit (1);
    }

    print_data (in);

    // 1. make the test_input.bin bigger so that we can actually test performance of encrpting raw data
    // 2. compare step 1 with the case if we encrypt huffman tree only in terms of performance
    // 3. understand what the compressed data look like....figure out what is the weakness if the raw compressed data is not encrypted.
    // 4. Revise the memory allocation code.


    // end of decription
#endif
    free(fprs[0]);
    free(fprs[1]);
    delete_test_key(ctx,key[0]);
    delete_test_key(ctx,key[1]);
    gpgme_data_release (in);
    gpgme_data_release (out);
    gpgme_release (ctx);
    // end of encryption


    // Write to output file
    writefile(argv[2], dest, destsize_bytes);

    // Free allocated memory
    free(src);
    free(dest);

    return 0;
}

// TODO: change get_frq_count so that it return max as well.
int quantize(float* src, size_t size, float error, int* quantized_src, float* min, float* max) {

    int frq_count = get_frq_count(src, size, error, min, max);
    
    float interval = *max - *min;
    float bucket_size = interval / frq_count;

    // Assign each data point to a bucket
    for (size_t i = 0; i < size; i++) {
        int bucket = static_cast<int>((src[i] - *min) / bucket_size);
        if (bucket >= frq_count) {
            bucket = frq_count - 1; // Handle edge case
        }

        // Store the bucket index
        quantized_src[i] = bucket;
    }

    return frq_count;
}

#if 0
int quantize(float* src, size_t size, float error, int* quantized_src, float* min) {
    // Find min and max values in data
    *min = src[0];
    float max = src[0];
    for (size_t i = 1; i < size; i++) {
        if (src[i] < *min)
            *min = src[i];
        if (src[i] > max)
            max = src[i];
    }

    // Calculate the interval and bucket width
    float interval = max - *min;
    if (interval == 0) {
        // All values are the same
        for (size_t i = 0; i < size; i++) {
            quantized_src[i] = 0;
        }
        return 1; // Only one bucket
    }
    int frq_count = static_cast<int>(interval / (error * 2) + 1); //frq_count as in number of buckets
    float bucket_size = interval / frq_count;


    // Assign each data point to a bucket
    for (size_t i = 0; i < size; i++) {
        int bucket = static_cast<int>((src[i] - *min) / bucket_size);
        if (bucket >= frq_count) {
            bucket = frq_count - 1; // Handle edge case
        }

        // Store the bucket index
        quantized_src[i] = bucket;
    }

    return frq_count;
}
#endif

void record_frequencies(const int* quantized_src, size_t fsize, size_t num_buckets, int* frqs) {
    // Count frequencies
    for (size_t i = 0; i < fsize; i++) {
        int bucket = quantized_src[i];
        if (bucket >= 0 && static_cast<size_t>(bucket) < num_buckets) {
            frqs[bucket]++; // Increment the frequency for the bucket
        }
    }
}

void make_queue(std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt>& phtree, int* frqs, size_t num_buckets) {
    // Iterate through all buckets
    for (size_t i = 0; i < num_buckets; ++i) {
        if (frqs[i] != 0) {
            // Create a new node with the bucket index and its frequency
            Node* new_node = new Node(static_cast<int>(i), frqs[i]);
            // Push the pointer to the node into the priority queue
            phtree.push(new_node);
        }
    }
}

Node* build_tree(std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt>& tree) {
    // Continue until there is only one node left in the priority queue
    while (tree.size() > 1) {
        // Extract the two nodes with the smallest frequencies
        Node* left = tree.top();
        tree.pop();
        Node* right = tree.top();
        tree.pop();

        // Create a new parent node with a combined frequency
        Node* parent = new Node(-1, left->freq + right->freq);
        parent->left = left;
        parent->right = right;

        // Push the new parent node back into the priority queue
        tree.push(parent);
    }

    // The last remaining node is the root of the Huffman tree
    return tree.top();
}

void assign_encode_helper(Node* node, unsigned int encode, int length) {
    if (NULL == node) {
        return;
    }

    // If the node is a leaf node (no left or right children)
    if (NULL == node->left && NULL == node->right) {
        // Assign the encoding and its length to the leaf node
        node->encode = encode;
        node->encode_length = length;
    } else {
        // Traverse the left subtree; append 0 to the encode
        assign_encode_helper(node->left, (encode << 1), length + 1);

        // Traverse the right subtree; append 1 to the encode
        assign_encode_helper(node->right, (encode << 1) | 1, length + 1);
    }
}

void assign_encode(Node* root) {
    if (NULL != root) {
        assign_encode_helper(root, 0, 0);
    }
}

void store_encodings_helper(Node* root, int* encodings, size_t frq_count) {
    if (!root)
        return;

    if (root->index >= 0 && static_cast<size_t>(root->index) < frq_count) {
        encodings[root->index * 2] = root->encode;
        encodings[root->index * 2 + 1] = root->encode_length;
    }

    if (root->left) {
        store_encodings_helper(root->left, encodings, frq_count);
    }

    if (root->right) {
        store_encodings_helper(root->right, encodings, frq_count);
    }
}

void store_encodings(Node* root, int* encodings, size_t frq_count) {
    // Initialize encodings
    // During decompression, you read 'encodings' from the buffer.
    // But  after that, you zero
    for (size_t i = 0; i < frq_count; ++i) {
        encodings[i * 2] = 0; // Initialize encode
        encodings[i * 2 + 1] = 0; // Initialize encode length
    }
    store_encodings_helper(root, encodings, frq_count);
}

void compress_quantization_levels(const int* src, char* dest, size_t num_elements, int* destsize, int* encodings) {
    int dest_index = 0;           // Tracks current position in the destination buffer
    unsigned int buffer = 0;      // Stores bits before writing to the destination buffer
    size_t src_index = 0;         // Tracks current position in the source buffer
    int bits_in_buffer = 0;       // Counts bits currently in the buffer
    *destsize = 0;                 // Size of the destination in bits

    // Iterate over each element in the source array
    while (src_index < num_elements) {
        int bucket_index = src[src_index];
        int encode = encodings[bucket_index * 2];       // Get Huffman encoding for the value
        int length = encodings[bucket_index * 2 + 1];   // Get the length of the encoding

        // Add the current value's encoding to the bit buffer
        buffer = (buffer << length) | encode;
        bits_in_buffer += length;

        // Write to the destination buffer byte by byte until there are less than 8 bits left
        while (bits_in_buffer >= 8) {
            unsigned char byte = (buffer >> (bits_in_buffer - 8)) & 0xFF;
            dest[dest_index++] = byte;    // Add the byte to the destination buffer
            bits_in_buffer -= 8;          // Update the number of bits in the buffer
            *destsize += 8;
        }
        src_index++;
    }

    // Handle any remaining bits in the buffer that don't fit into a full byte
    if (bits_in_buffer > 0) {
        unsigned char byte = buffer << (8 - bits_in_buffer);
        dest[dest_index++] = byte;
        *destsize += bits_in_buffer;
    }
    // No need to null-terminate binary data
}

void reverse_quantize(const int* quantized_data, size_t size, float min, float bucket_size, float* data) {
    for (size_t i = 0; i < size; i++) {
        // Calculate the approximate original value by finding the midpoint of the bucket
        data[i] = min + quantized_data[i] * bucket_size + (bucket_size / 2.0f);
    }
}

double calc_speed(long original_size, double compression_time) {
    if (compression_time == 0) {
        return 0;
    }
    return original_size / (compression_time * 1024 * 1024);
}

void free_tree(Node* root) {
    if (!root)
        return;

    // Recursively free left and right subtrees
    if (root->left) {
        free_tree(root->left);
    }

    if (root->right) {
        free_tree(root->right);
    }

    // Free the current node
    delete root; // Use delete instead of free for C++ objects
}

// int compress(float* src, char* dest, size_t num_elements, int* destsize, float error) {
//     // Allocate memory for quantized_src
//     int* quantized_src = (int*)malloc(num_elements * sizeof(int));
//     if (NULL == quantized_src) {
//         printf("Error allocating memory for quantized_src\n");
//         return 1;
//     }

//     // Quantization
//     float min = 0.0f;
//     float bucket_size = 0.0f;
//     int frq_count = quantize(src, num_elements, error, quantized_src, &min);
//     bucket_size = error * 2; // As per your original code

//     // Allocate and initialize frequency array
//     int* frqs = (int*)malloc(frq_count * sizeof(int));
//     if (NULL == frqs) {
//         printf("Error allocating memory for frqs\n");
//         free(quantized_src);
//         return 1;
//     }
//     memset(frqs, 0, frq_count * sizeof(int));

//     // Record frequencies
//     record_frequencies(quantized_src, num_elements, frq_count, frqs);

//     // Build Huffman Tree
//     std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt> tree;
//     make_queue(tree, frqs, frq_count);
//     Node* root = build_tree(tree);
//     assign_encode(root);

//     // Allocate encodings array as a contiguous block
//     int* encodings = (int*)malloc(frq_count * 2 * sizeof(int));
//     if (NULL == encodings) {
//         printf("Error allocating memory for encodings\n");
//         free(quantized_src);
//         free(frqs);
//         free_tree(root);
//         return 1;
//     }

//     // Store encodings
//     store_encodings(root, encodings, frq_count);

//     // Perform compression
//     compress_quantization_levels(quantized_src, dest, num_elements, destsize, encodings);

//     // Free allocated memory except frqs and encodings (needed for serialization)
//     free(quantized_src);
//     free_tree(root);

//     // Serialize metadata into the dest buffer
//     size_t header_size = sizeof(size_t) + sizeof(float) + sizeof(float) + sizeof(int) +
//                          frq_count * sizeof(int) + frq_count * 2 * sizeof(int);
//     size_t total_size = header_size + ((*destsize + 7) / 8);

//     // Shift existing data in dest to make space for the header
//     memmove(dest + header_size, dest, (*destsize + 7) / 8);

//     // Serialize header into dest buffer
//     size_t offset = 0;
//     memcpy(dest + offset, &num_elements, sizeof(size_t));
//     offset += sizeof(size_t);
//     memcpy(dest + offset, &min, sizeof(float));
//     offset += sizeof(float);
//     memcpy(dest + offset, &bucket_size, sizeof(float));
//     offset += sizeof(float);
//     memcpy(dest + offset, &frq_count, sizeof(int));
//     offset += sizeof(int);
//     memcpy(dest + offset, frqs, frq_count * sizeof(int));
//     offset += frq_count * sizeof(int);
//     memcpy(dest + offset, encodings, frq_count * 2 * sizeof(int));
//     offset += frq_count * 2 * sizeof(int);

   
//     printf("numelems: %d\n", num_elements);
//     printf("min: %f\n", min);
//     printf("bucket_size: %f\n", bucket_size);
//     printf("frq_count: %d\n", frq_count);
//     printf("offset: %d\n", offset);
//     // Update destsize to include header
//     *destsize = static_cast<int>(total_size);

//     // Free frqs and encodings as they are now serialized into dest
//     free(frqs);
//     free(encodings);

//     return 0;
// }


int compress(float  *src,
    char   *dest,
    size_t  num_elements,
    int    *destsize_bytes,
    float   error)
{
    /* ---------- quantisation ---------- */
    std::vector<int> qsrc(num_elements);
    float min = 0.f, max = 0.f, bucket_sz = 0.f;

    // 6‑argument quantize: returns frq_count and fills min & max
    const int frq_count =
    quantize(src, num_elements, error, qsrc.data(), &min, &max);

    bucket_sz = error * 2.f;              // bucket width = 2·error

    /* ---------- frequency table ---------- */
    std::vector<int> frqs(frq_count, 0);
    record_frequencies(qsrc.data(), num_elements, frq_count, frqs.data());

    /* ---------- Huffman tree ---------- */
    std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt> pq;
    make_queue(pq, frqs.data(), frq_count);
    Node *root = build_tree(pq);
    assign_encode(root);

    std::vector<int> enc(frq_count * 2);
    store_encodings(root, enc.data(), frq_count);

    /* ---------------------------------------------------------
    HEADER BYTES  (scalars + frqs[] only)
    --------------------------------------------------------- */
    //changed variable name header_bytes
    //Gary
    const size_t header_bytes = sizeof(size_t)          /* num_elements */
                        + 2 * sizeof(float)       /* min, bucket  */
                        + sizeof(int)             /* frq_count    */
                        + frq_count * sizeof(int);/* frqs[]       */

    /* ---------- write bit‑stream right after header ---------- */
    int bit_cnt = 0;
    compress_quantization_levels(qsrc.data(),
                            dest + header_bytes,
                            num_elements,
                            &bit_cnt,
                            enc.data());
    const size_t payload_bytes = (bit_cnt + 7) / 8;      // ceil(bits/8)

    /* ---------- serialize header ---------- */
    char *cur = dest;
    auto copy_scalar = [&](auto v) {
    cur = std::copy_n(reinterpret_cast<const char*>(&v), sizeof(v), cur);
    };
    copy_scalar(num_elements);
    copy_scalar(min);
    copy_scalar(bucket_sz);
    copy_scalar(frq_count);

    cur = std::copy_n(reinterpret_cast<const char*>(frqs.data()),
                frq_count * sizeof(int),
                cur);

    /* (We no longer store enc[] in the header; decompressor will rebuild) */

    free_tree(root);
    *destsize_bytes = static_cast<int>(header_bytes + payload_bytes);
    return 0;
}

void writefile(const char* fname, const char* buff, size_t size) {
    FILE* output = fopen(fname, "wb");
    if (NULL == output) {
        perror("Error opening output file");
        return;
    }

    size_t written = fwrite(buff, 1, size, output);
    if (written != size) {
        printf("Error writing to output file\n");
    }

    fclose(output);
}


