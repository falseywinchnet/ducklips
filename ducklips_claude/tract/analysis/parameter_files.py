"""Block files of the engine's parameter tree (tract/parameter_tree.hpp):
read as [dotted path, value] in file order, written back grouped by block."""


def split_words(text):
    words = []
    for line in text.splitlines():
        line = line.split("#", 1)[0].replace("{", " { ").replace("}", " } ")
        words.extend(line.split())
    return words


def read_tree(path):
    """The leaves of a block file as [dotted path, value], in file order."""
    with open(path) as handle:
        words = split_words(handle.read())
    entries = []
    open_blocks = []
    at = 0
    while at < len(words):
        word = words[at]
        if word == "}":
            open_blocks.pop()
            at += 1
        elif at + 1 < len(words) and words[at + 1] == "{":
            open_blocks.append(word)
            at += 2
        else:
            entries.append([".".join(open_blocks + [word]), float(words[at + 1])])
            at += 2
    return entries


def write_tree(path, entries, header=""):
    lines = [header.rstrip("\n")] if header else []
    open_blocks = []
    for name, value in entries:
        parts = name.split(".")
        shared = 0
        while shared < len(open_blocks) and shared + 1 < len(parts) and open_blocks[shared] == parts[shared]:
            shared += 1
        while len(open_blocks) > shared:
            open_blocks.pop()
            lines.append("  " * len(open_blocks) + "}")
        while len(open_blocks) + 1 < len(parts):
            lines.append("  " * len(open_blocks) + parts[len(open_blocks)] + " {")
            open_blocks.append(parts[len(open_blocks)])
        lines.append("  " * len(open_blocks) + "%s %.10g" % (parts[-1], value))
    while open_blocks:
        open_blocks.pop()
        lines.append("  " * len(open_blocks) + "}")
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")
