def to_camel(string: str) -> str:
    words = [word for word in string.split('_')]
    return words[0] + ''.join(w.capitalize() for w in words[1:])
