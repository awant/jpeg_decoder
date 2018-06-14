#pragma once

#include <vector>
#include <iostream>
#include <assert.h>

template <class T>
struct Node {
    T value = 0;
    std::unique_ptr<Node<T>> left = nullptr;
    std::unique_ptr<Node<T>> right = nullptr;
    Node<T>* parent = nullptr;
};

template <class T>
class HuffmanTree {
public:
    const int MAX_HUFFMAN_COEFF_SIZE = 16;

    class Iterator {
    public:
        explicit Iterator(Node<T>* node)
                : current_node_(node) {}

        void LeftStep() {
            assert(current_node_ != nullptr);
            current_node_ = current_node_->left.get();
        }

        void RightStep() {
            assert(current_node_ != nullptr);
            current_node_ = current_node_->right.get();
        }

        bool operator != (const Iterator& it) {
            return current_node_ != it.current_node_;
        }

        bool Last() {
            assert(current_node_ != nullptr);
            return (current_node_->left.get() == nullptr) and
                   (current_node_->right.get() == nullptr);
        }

        T& operator*() {
            return current_node_->value;
        }

    private:
        Node<T>* current_node_ = nullptr;
    };

    HuffmanTree(): tree_(std::make_unique<Node<T>>()) {}

    HuffmanTree(HuffmanTree &&htree) noexcept {
        tree_ = std::move(htree.tree_);
    }

    HuffmanTree(const std::vector<int>& codes_counter,
                const std::vector<T>& codes);

    T Decode(const std::vector<char>& sequence, bool* is_decoded);

    Iterator Begin() {
        return Iterator(tree_.get());
    }


    void Print();
private:
    void Print(Node<T>* node, const std::string& cur_path);

    std::unique_ptr<Node<T>> tree_ = nullptr;

    Node<T>* current_node_ = nullptr;

    void SetNewNode(const T& val, int len);
};

template<class T>
HuffmanTree<T>::HuffmanTree(const std::vector<int>& codes_counter,
            const std::vector<T>& codes) {
    tree_.reset(new Node<T>);

    if (codes_counter.size() != MAX_HUFFMAN_COEFF_SIZE) {
        throw std::runtime_error("Huffman codes counter != " + std::to_string(MAX_HUFFMAN_COEFF_SIZE));
    }
    int current_length = 0;
    current_node_ = tree_.get();
    int code_idx = 0;
    for (int code_length = 1; code_length <= 16; ++code_length) {
        int code_counter = codes_counter[code_length-1];
        for (int i = 0; i < code_counter; ++i) {
            SetNewNode(codes[code_idx++], code_length - current_length);
            current_length = code_length - 1;
        }
    }
}

template<class T>
T HuffmanTree<T>::Decode(const std::vector<char>& sequence, bool* is_decoded) {
    current_node_ = tree_.get();
    for (const auto& val : sequence) {
        if ((val == 0) and (current_node_->left)) {
            current_node_ = current_node_->left.get();
        } else if ((val == 1) and (current_node_->right)) {
            current_node_ = current_node_->right.get();
        } else {
            throw std::runtime_error("Huffman decoding: wrong direction");
        }
    }
    *is_decoded = (!current_node_->left) and (!current_node_->right);
    return current_node_->value;
}

template<class T>
void HuffmanTree<T>::SetNewNode(const T& val, int len) {
    if (current_node_->left) {
        if (current_node_->right) {
            current_node_ = current_node_->parent;
            SetNewNode(val, len+1);
            return;
        }
        auto* node = new Node<T>;
        node->parent = current_node_;
        current_node_->right.reset(node);
        if (len == 1) {
            node->value = val;
            return;
        }
        current_node_ = node;
        SetNewNode(val, len-1);
    } else {
        auto* node = new Node<T>;
        node->parent = current_node_;
        current_node_->left.reset(node);
        if (len == 1) {
            node->value = val;
            return;
        }
        current_node_ = node;
        SetNewNode(val, len-1);
    }
}
