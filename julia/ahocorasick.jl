module AhoCorasick

import Base.show

type Node
    depth :: Int
    is_terminal :: Bool
    children :: Dict{Char, Node}
    fail :: Union(Nothing, Node)
end

function Base.show(io :: IO, node :: Node)
    write(io, "Node(depth=$(node.depth))")
end 
 
type Automaton
    root :: Node
    case_sensitive :: Bool
end

function new_node(depth :: Int = 0, is_terminal :: Bool = false)
    Node(depth, is_terminal, Dict{Char,Node}(), nothing)
end

function new(case_sensitive :: Bool = true)
    Automaton(new_node(), case_sensitive)
end

function add(node :: Node, word :: UTF8String)
    if length (word) == 0
        node.is_terminal = true
    else
        c = word[1]
        if !haskey(node.children, c)
            node.children[c] = new_node(node.depth + 1)
        end
        add(node.children[c], word[2:end])
    end
end

function search (node :: Node, c :: Char)
    if haskey (node.children, c)
        node.children[c]
    elseif node.depth == 0
        node
    elseif node.fail == nothing
        node
    else
        search(node.fail, c)
    end
end

function add_fail_transition!(node :: Node)
    for (c,child) in node.children
        child.fail = search(node.fail, c)
    end
    for child in values(node.children)
        add_fail_transition!(child)
    end
end

function build(ac :: Automaton)
    ac.root.fail = ac.root
    for child in values (ac.root.children)
        child.fail = ac.root
    end
    for child in values (ac.root.children)
        add_fail_transition!(child)
    end
end
 
function search (ac :: Automaton, text :: UTF8String)
    if !ac.case_sensitive
        text = lowercase (text)
    end
    i = 0
    node = ac.root
    matches = Range1{Int}[]
    while length (text) > 0
        i += 1
        node = search (node, text[1])
        if node.is_terminal
            loc = (1+i-node.depth):i
            push!(matches, loc)
        end
        text = text[2:end]
    end
    matches
end

function search (ac :: Automaton, text :: String)
    search(ac, utf8(text))
end

function add (ac :: Automaton, word :: UTF8String)
    if !ac.case_sensitive
        text = lowercase(word)
    end
    add (ac.root, word)
end

function add (ac :: Automaton, word :: String)
    add (ac, utf8(word))
end

end

ac = AhoCorasick.new(false)
AhoCorasick.add(ac, "alpha")
AhoCorasick.add(ac, "beta")
AhoCorasick.build(ac)
text = "the Alpha glycoprotein beta and alpha time"
matches = AhoCorasick.search(ac, text)
map(x->text[x], matches)
