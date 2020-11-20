$('a').click(function (e) {
    var x = window.pageXOffset,
        y = window.pageYOffset;
    $(window).one('scroll', function () {
        window.scrollTo(x, y);
    })
});