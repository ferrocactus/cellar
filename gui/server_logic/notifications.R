notifications <- function(input, output, session) {
    notification_list = reactiveVal()
    output$notifications = renderMenu({
        validate(need(notification_list(), message = FALSE))
        dropdownMenu(type = 'notifications',
                     badgeStatus = 'warning', .list = notification_list())
    })

    return(list(
        push = function(message) {
            pf = parent.env(environment())
            pf$notification_list(c(pf$notification_list(), list(message)))
        },
        pop = function() {
            pf = parent.env(environment())
            pf$notification_list(
                notification_list()[-length(pf$notification_list())])
        }
  ))
}
