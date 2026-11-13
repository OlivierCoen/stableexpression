import nltk
from nltk.corpus import wordnet
from functools import lru_cache

nltk.download("punkt_tab")
nltk.download("averaged_perceptron_tagger_eng")
nltk.download("wordnet")

lemmatizer = nltk.WordNetLemmatizer()
stemmer = nltk.PorterStemmer()


def get_wordnet_pos(token: str) -> str:
    """Map POS tag to first character lemmatize() accepts."""
    tag = nltk.pos_tag([token])[0][1][0].upper()
    tag_dict = {
        "J": wordnet.ADJ,
        "N": wordnet.NOUN,
        "V": wordnet.VERB,
        "R": wordnet.ADV,
    }
    return tag_dict.get(tag, wordnet.NOUN)  # Default to NOUN if not found


@lru_cache
def get_stemmed_tokens(sentence: str) -> list[str]:
    """Tokenize a sentence into its constituent words, and then stem each word.

    Parameters
    ----------
    sentence : str
        The sentence to be tokenized and stemmed

    Returns
    -------
    tokens : List[str]
        The list of stemmed tokens

    """
    tokens = nltk.word_tokenize(sentence)
    return [stemmer.stem(token) for token in tokens]


@lru_cache
def get_lemmed_tokens(sentence: str) -> list[str]:
    """Tokenize a sentence into its constituent words, and then lemmatize each word.

    Parameters
    ----------
    sentence : str
        The sentence to be tokenized and lemmatized

    Returns
    -------
    tokens : List[str]
        The list of lemmatized tokens

    """
    tokens = nltk.word_tokenize(sentence)
    return [lemmatizer.lemmatize(token, get_wordnet_pos(token)) for token in tokens]


@lru_cache
def get_synonyms(word: str) -> set[str]:
    """Get all synonyms of a word from the wordnet database.

    Parameters
    ----------
    word : str
        The word for which to get synonyms

    Returns
    -------
    synonyms : set
        A set of all synonyms of the word

    """
    synonyms = []
    for syn in wordnet.synsets(word):
        for lemma in syn.lemmas():
            synonyms.append(lemma.name())  # Get the name of each lemma (synonym)
    return set(synonyms)  # Return as a set to avoid duplicates


@lru_cache
def get_all_candidate_target_words(sentence: str) -> list[str]:
    """Get all candidate target words from a sentence by stemming and lemmatizing the tokens and getting synonyms from the wordnet database.

    Parameters
    ----------
    sentence : str
        The sentence from which to get candidate target words

    Returns
    -------
    candidates : list
        A list of all candidate target words

    """
    candidates = []
    lemmatized_tokens = get_lemmed_tokens(sentence)
    stemmed_tokens = get_stemmed_tokens(sentence)
    tokens = list(set(lemmatized_tokens + stemmed_tokens))
    for token in tokens:
        candidates += get_synonyms(token)
    return candidates


def keywords_in_fields(fields: list[str], keywords: list[str]) -> list[str]:
    """Check if any of the keywords are in any of the fields.

    Parameters
    ----------
    fields : list[str]
        The fields to be searched
    keywords : list[str]
        The keywords to be searched for

    Returns
    -------
    list[str]
        A list of all keywords that are found in any of the fields

    """
    matched = set()
    keyword_stems = {kw: frozenset([kw, *get_stemmed_tokens(kw)]) for kw in keywords}
    for field in fields:
        candidates = get_all_candidate_target_words(field)  # cached, so once per distinct field
        for kw, stems in keyword_stems.items():
            if any(s in field or any(s in c for c in candidates) for s in stems):
                matched.add(kw)
    return list(matched)
