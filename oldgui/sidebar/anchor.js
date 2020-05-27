window.onload = function() {
    var anchors = document.getElementsByTagName("a");

    for (var i = 0; i < anchors.length; i++) {
        if (anchors[i].href.endsWith("#")) {
            anchors[i].href = anchors[i].href + "!"
        }
    }
}
